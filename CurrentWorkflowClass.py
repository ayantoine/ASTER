# coding: utf-8
"""Python3.6"""

import os
import sys
import re
import copy

import inspect

from Bio.Blast import NCBIXML
from Bio import SearchIO
from Bio import SeqIO

########################################################################
#CONSTANT
FEATURE_TRANSCRIPT="transcript"
FEATURE_EXON="exon"
FEATURE_GENE="gene"
STRAND_TO_INT={"+":1,"-":-1}

EXON_BASENAME="Ex"
EXON_SEPARATOR="."
EXON_INTERNALSEPARATOIR="-"

BSTDOUT=True

DEBUG_BOOL=False
DEBUG_READID="ch117_read1370_template_pass_BYK_CB_ONT_1_FAF04998_A"

#Exonerate Parsing
REVERSE_COMP={"A":"T","C":"G","G":"C","T":"A",
                "a":"t","c":"g","g":"c","t":"a",
                "N":"N","n":"n",
                ">":"<","<":">",".":".","-":"-"," ":" "
                }
RAWSCORE_THRESHOLD=150

REGROUP_READEXON_VALUE_BOOL=True
REGROUP_READEXON_VALUE=6
REGROUP_READEXON_THRESHOLD=0.1 #if not REGROUP_READEXON_VALUE_BOOL, use REGROUP_READEXON_THRESHOLD

REGROUP_LASTREADEXON_THRESHOLD=0.1
REGROUP_FIRSTREADEXON_THRESHOLD=0.1

REGROUP_GLOBALVECTOR_VALUE=6

GFF_COLUMN=["seqname","source","feature","start","end","score","strand","frame","attribute"]
PAF_COLUMN=["ReadId","ReadSize","ReadStart","ReadStop","Strand","GeneId","GeneSize","GeneStart","GeneStop","NbrNotIndel","AlignSize","Quality"]
########################################################################
#Function
def Fasta2Dict(sFile):
    dResults={}
    sSeqName=""
    sSeqContent=""
    for sLine in open(sFile):
        sStripedLine=sLine.strip()
        if sStripedLine[0]==">":
            if sSeqName!="":
                dResults[sSeqName]=sSeqContent
            sSeqName=sStripedLine[1:]
            sSeqContent=""
        else:
            sSeqContent+=sStripedLine
    if sSeqName!="":
        dResults[sSeqName]=sSeqContent
    return dResults

########################################################################
#CLASS
class PseudoBlastMatch:
    
    def __init__(self,sGeneId,iGeneSize):
        self.query=sGeneId
        self.query_letters=iGeneSize
        self.alignments=[]

class PseudoBlastAlignments:

    def __init__(self,sReadId,iReadSize):
        self.title=sReadId
        self.length=iReadSize
        self.hsps=[]
        
class PseudoBlastHSP:
    
    def __init__(self,sGenePart,iGenePartStart,iGenePartEnd,sReadPart,iReadPartStart,iReadPartEnd):
        self.query=sGenePart
        self.query_start=iGenePartStart
        self.query_end=iGenePartEnd
        self.sbjct=sReadPart
        self.sbjct_start=iReadPartStart
        self.sbjct_end=iReadPartEnd
        
        if iGenePartEnd>iGenePartStart:
            iGeneStrand=1
        else:
            iGeneStrand=-1
        if iReadPartEnd>iReadPartStart:
            iReadStrand=1
        else:
            iReadStrand=-1
        self.frame=(iGeneStrand,iReadStrand)
        self.strand=(None,None)


class ToolContent:
    def __init__(self):
        self.gene_strand=None
        self.tool=None
        
    def set_tool(self,sString):
        self.tool=sString
    
    def set_queryId(self, sString):
        self.query_id=sString
    
    def get_queryId(self):
        return self.query_id
    
    def set_querySize(self, sString):
        self.query_size=sString
        
    def set_emptyAlignList(self):
        self.align_list=[]
        
    def add_align(self,oAlign):
        if oAlign.get_read_strand() is not None and oAlign.get_gene_strand() is not None:
            self.align_list.append(oAlign)
    
    def set_alignList(self,tObjectList):
        self.align_list=list(tObjectList)
        
    def get_alignList(self):
        return self.align_list
        
    def get_strand(self):
        return self.strand
        
    def set_strand(self,iInt):
        self.strand=iInt
        
    def compute_self_cover_alignments(self):
        for oAlignContent in self.get_alignList():
            oAlignContent.compute_self_cover()
            
    def compute_ref_cover_alignments(self):
        for oAlignContent in self.get_alignList():
            oAlignContent.compute_ref_cover(self.get_querySize())
            
    def define_position_alignments(self):
        for oAlignContent in self.get_alignList():
            oAlignContent.define_position()
        self.compute_ref_cover_alignments()
        
    def describe_AllHsp(self):
        for oAlignContent in self.get_alignList():
            print("=======================")
            print(oAlignContent.get_id())
            print("---------")
            for oHsp in oAlignContent.get_hspList():
                print(oHsp.sbjct_start,oHsp.sbjct_end,len(oHsp.sbjct.replace("-","")))
                print(oHsp.query_start,oHsp.query_end,len(oHsp.query.replace("-","")))
                print("---------")
    
    def describe_alignment(self,oTranscriptContent=None,bPrint=True):
        sContent=""
        for oAlignContent in self.get_alignList():
            sContent+="{}\t".format(self.tool)
            sContent+=oAlignContent.describe_self(oTrContent=oTranscriptContent,bStdout=bPrint)
        return sContent
        
    def describe_specificAlignment(self,sId,oTranscriptContent=None,bPrint=True):
        sContent=""
        sContent+="{}\t".format(self.tool)
        bWrited=False
        for oAlignContent in self.get_alignList():
            if oAlignContent.get_id()==sId:
                sContent+=oAlignContent.describe_self(oTrContent=oTranscriptContent,bStdout=bPrint)
                bWrited=True
                break
        if not bWrited:
            sContent+=oAlignContent.describe_empty(bStdout=bPrint)
        return sContent
    
    def describe_alignment_header(self,oTranscriptContent=None,bPrint=True):
        sContent=""
        for oAlignContent in self.get_alignList():
            sContent+=oAlignContent.describe_header(oTrContent=oTranscriptContent,bStdout=bPrint)
            break
        return sContent
    
    def get_specificAlignment_alignContent_size(self,sId):
        for oAlignContent in self.get_alignList():
            if oAlignContent.get_id()==sId:
                return oAlignContent.get_alignContent_selfCover()
        return 0.0
        
    def assign_exonCovering(self,oTranscriptContent):
        for oAlignContent in self.get_alignList():
            oAlignContent.assign_exonCovering(oTranscriptContent)
            
    def purge_alignList(self,tTargetIdToKeep):
        tNewList=[]
        for oAlignContent in self.get_alignList():
            if oAlignContent.get_id() in tTargetIdToKeep:
                tNewList.append(oAlignContent)
        self.set_alignList(tNewList)
        
    def get_querySize(self):
        return self.query_size
        
    def get_alignId(self):
        tList=[]
        for oAlignContent in self.get_alignList():
            tList.append(oAlignContent.get_id())
        return tList
    
    def export_data(self):
        sContent=""
        sGeneId=self.get_queryId()
        iGeneSize=self.get_querySize()

        for oAlignContent in self.get_alignList():
            sReadId=oAlignContent.get_id()
            iReadSize=oAlignContent.get_read_size()

            for oHsp in oAlignContent.get_hspList():
                iReadStart=min(oHsp.sbjct_start,oHsp.sbjct_end)
                iReadStop=max(oHsp.sbjct_start,oHsp.sbjct_end)
                iGeneStart=min(oHsp.query_start,oHsp.query_end)
                iGeneStop=max(oHsp.query_start,oHsp.query_end)
                iReadStrand=oHsp.frame[1]
                iGeneStrand=oHsp.frame[0]
                iAlignSize=len(oHsp.sbjct)
                iNumberExactMatch=0
                iNumberIndel=0
                iNumberNotIndel=0
                for iIndex in range(iAlignSize):
                    if oHsp.sbjct[iIndex]==oHsp.query[iIndex]:
                        iNumberExactMatch+=1
                    if oHsp.sbjct[iIndex]=="-" or oHsp.query[iIndex]=="-":
                        iNumberIndel+=1
                    else:
                        iNumberNotIndel+=1
                if iReadStrand*iGeneStrand==1:
                    sReadStrand="+"
                elif iReadStrand*iGeneStrand==-1:
                    sReadStrand="-"
                else:
                    exit("FATAL 235")
                tFormatArgs=[sReadId,iReadSize,iReadStart,iReadStop,sReadStrand,
                            sGeneId,iGeneSize,iGeneStart,iGeneStop,
                            iNumberNotIndel,iAlignSize,
                            255]
                sContent+="\t".join(["{}"]*len(tFormatArgs)).format(*tFormatArgs)
                sContent+="\n"
        return sContent
                    
    """
    Two Hsp coord that overlap on the Read must overlap on the Gene
    Two Hsp coord that overlap on the Gene must overlap on the Read 
    """
    def check_divergence(self):
        for oAlignContent in self.get_alignList():
            tHspList=oAlignContent.get_hspList()
            while True:
                tRemoveHsp=[]
                for iCurrentHspIndex in range(len(tHspList)-1):
                    oCurrentHsp=tHspList[iCurrentHspIndex]
                    iCurrentReadStart=min(oCurrentHsp.sbjct_start,oCurrentHsp.sbjct_end)
                    iCurrentReadStop=max(oCurrentHsp.sbjct_start,oCurrentHsp.sbjct_end)
                    iCurrentGeneStart=min(oCurrentHsp.query_start,oCurrentHsp.query_end)
                    iCurrentGeneStop=max(oCurrentHsp.query_start,oCurrentHsp.query_end)
                    
                    iCurrentSize=len(oCurrentHsp.sbjct)
                    iCurrentError=0
                    for iIndex in range(len(oCurrentHsp.sbjct)):
                        if oCurrentHsp.sbjct[iIndex]!=oCurrentHsp.query[iIndex]:
                            iCurrentError+=1
                    fCurrentRatioError=float(iCurrentError)/iCurrentSize
                    
                    for iAnotherHspIndex in range(iCurrentHspIndex+1,len(tHspList)):
                        oAnotherHsp=tHspList[iAnotherHspIndex]
                        iAnotherReadStart=min(oAnotherHsp.sbjct_start,oAnotherHsp.sbjct_end)
                        iAnotherReadStop=max(oAnotherHsp.sbjct_start,oAnotherHsp.sbjct_end)
                        iAnotherGeneStart=min(oAnotherHsp.query_start,oAnotherHsp.query_end)
                        iAnotherGeneStop=max(oAnotherHsp.query_start,oAnotherHsp.query_end)
                        
                        iAnotherSize=len(oAnotherHsp.sbjct)
                        iAnotherError=0
                        for iIndex in range(len(oAnotherHsp.sbjct)):
                            if oAnotherHsp.sbjct[iIndex]!=oAnotherHsp.query[iIndex]:
                                iAnotherError+=1
                        fAnotherRatioError=float(iAnotherError)/iAnotherSize                        
                        
                        bCurrentHaveBestAlignment=True
                        if fAnotherRatioError<fCurrentRatioError:
                            bCurrentHaveBestAlignment=False
                        
                        if iAnotherReadStart<=iCurrentReadStop and iAnotherReadStop>=iCurrentReadStart:
                            #Check if coord on Read are overlaped
                            if iAnotherReadStart<=iCurrentReadStart and iAnotherReadStop>=iCurrentReadStop:
                                #Current is include into Another
                                if iAnotherGeneStart<=iCurrentGeneStart and iAnotherGeneStop>=iCurrentGeneStop:
                                    #Same thing in Gene. Remove Current
                                    tRemoveHsp.append(oCurrentHsp)
                                    break
                                else:
                                    #Not the same thing in Gene. Keep only the best alignment
                                    if bCurrentHaveBestAlignment:
                                        #print("Current_sbjct",oCurrentHsp.sbjct_start,oCurrentHsp.sbjct_end,len(oCurrentHsp.sbjct.replace("-","")))
                                        #print("Current_query",oCurrentHsp.query_start,oCurrentHsp.query_end,len(oCurrentHsp.query.replace("-","")))
                                        #print(iCurrentSize,iCurrentError,fCurrentRatioError)
                                        #print("Another_sbjct",oAnotherHsp.sbjct_start,oAnotherHsp.sbjct_end,len(oAnotherHsp.sbjct.replace("-","")))
                                        #print("Another_query",oAnotherHsp.query_start,oAnotherHsp.query_end,len(oAnotherHsp.query.replace("-","")))
                                        #print(iAnotherSize,iAnotherError,fAnotherRatioError)
                                        tRemoveHsp.append(oAnotherHsp)
                                    else:
                                        tRemoveHsp.append(oCurrentHsp)
                                    break
                            elif iAnotherReadStart>=iCurrentReadStart and iAnotherReadStop<=iCurrentReadStop:
                                #Another is include into Current
                                if iAnotherGeneStart>=iCurrentGeneStart and iAnotherGeneStop<=iCurrentGeneStop:
                                    #Same thing in Gene. Remove Another
                                    tRemoveHsp.append(oAnotherHsp)
                                    break
                                else:
                                    #Not the same thing in Gene. Keep only the best alignment
                                    if bCurrentHaveBestAlignment:
                                        tRemoveHsp.append(oAnotherHsp)
                                    else:
                                        tRemoveHsp.append(oCurrentHsp)
                                    break
                            elif iAnotherReadStart<=iCurrentReadStart and iAnotherReadStop>=iCurrentReadStart:
                                #Another overlap the begining of Current
                                if iAnotherGeneStart<=iCurrentGeneStart and iAnotherGeneStop>=iCurrentGeneStart:
                                    #Same thing in Gene. Do nothing.
                                    pass
                                else:
                                    #Not the same thing in Gene. Kee only the best alignment
                                    if bCurrentHaveBestAlignment:
                                        tRemoveHsp.append(oAnotherHsp)
                                    else:
                                        tRemoveHsp.append(oCurrentHsp)
                                    break
                            elif iAnotherReadStart<=iCurrentReadStop and iAnotherReadStop>=iCurrentReadStop:
                                #Current overlap the ending of Another
                                if iAnotherGeneStart<=iCurrentGeneStop and iAnotherGeneStop>=iCurrentGeneStop:
                                    #Same thing in Gene. Do nothing.
                                    pass
                                else:
                                    #Not the same thing in Gene. Kee only the best alignment
                                    if bCurrentHaveBestAlignment:
                                        tRemoveHsp.append(oAnotherHsp)
                                    else:
                                        tRemoveHsp.append(oCurrentHsp)
                                    break
                        
                        if iAnotherGeneStart<=iCurrentGeneStop and iAnotherGeneStop>=iAnotherGeneStart:    
                            #Check if coord on Gene are overlaped
                            if iAnotherGeneStart<=iCurrentGeneStart and iAnotherGeneStop>=iCurrentGeneStop:
                                #Current is include into Another
                                if iAnotherReadStart<=iCurrentReadStart and iAnotherReadStop>=iCurrentReadStop:
                                    #Same thing in Read. Remove Current
                                    tRemoveHsp.append(oCurrentHsp)
                                    break
                                else:
                                    #Not the same thing in Gene. Kee only the best alignment
                                    if bCurrentHaveBestAlignment:
                                        tRemoveHsp.append(oAnotherHsp)
                                    else:
                                        tRemoveHsp.append(oCurrentHsp)
                                    break
                            elif iAnotherGeneStart>=iCurrentGeneStart and iAnotherGeneStop<=iCurrentGeneStop:
                                #Another is include into Current
                                if iAnotherReadStart>=iCurrentReadStart and iAnotherReadStop<=iCurrentReadStop:
                                    #Same thing in Read. Remove Current
                                    tRemoveHsp.append(oCurrentHsp)
                                    break
                                else:
                                    #Not the same thing in Gene. Kee only the best alignment
                                    if bCurrentHaveBestAlignment:
                                        tRemoveHsp.append(oAnotherHsp)
                                    else:
                                        tRemoveHsp.append(oCurrentHsp)
                                    break
                            elif iAnotherGeneStart<=iCurrentGeneStart and iAnotherGeneStop>=iCurrentGeneStart:
                                #Another overlap the begining of Current
                                if iAnotherReadStart<=iCurrentReadStart and iAnotherReadStop>=iCurrentReadStart:
                                    #Same thing in Read. Do nothing.
                                    pass
                                else:
                                    #Not the same thing in Gene. Kee only the best alignment
                                    if bCurrentHaveBestAlignment:
                                        tRemoveHsp.append(oAnotherHsp)
                                    else:
                                        tRemoveHsp.append(oCurrentHsp)
                                    break
                            elif iAnotherGeneStart>=iCurrentGeneStop and iAnotherGeneStop>=iCurrentGeneStop:
                                #Current overlap the ending of Another
                                if iAnotherReadStart<=iCurrentReadStop and iAnotherReadStop>=iCurrentReadStop:
                                    #Same thing in Gene. Do nothing.
                                    pass
                                else:
                                    #Not the same thing in Gene. Kee only the best alignment
                                    if bCurrentHaveBestAlignment:
                                        tRemoveHsp.append(oAnotherHsp)
                                    else:
                                        tRemoveHsp.append(oCurrentHsp)
                                    break

                    if len(tRemoveHsp)!=0:
                        break
                if len(tRemoveHsp)==0:
                    break
                else:
                    print("WARNING : {} reject {} Hsp".format(oAlignContent.get_id(),len(tRemoveHsp)))
                    #for oHsp in tRemoveHsp:
                        #print(oHsp.sbjct_start,oHsp.sbjct_end,len(oHsp.sbjct.replace("-","")))
                    tHspList=[X for X in tHspList if X not in tRemoveHsp]
                    if len(tHspList)<2:
                        break
            oAlignContent.set_hspList(tHspList)
                    
                    
    
    def merge_withContent(self,oToolContent):
        self.set_tool("merge")
        #Make fusion for AlignContent if exists in both ToolContent
        for oSelfAlignContent in self.get_alignList():
            sSelfAlignId=oSelfAlignContent.get_id()
            iSelfSize=oSelfAlignContent.get_alignContent_size()
            tNewHspList=[]
            tAllHspList=list(oSelfAlignContent.get_hspList())
            for oOtherAlignContent in oToolContent.get_alignList():
                sOtherAlignId=oOtherAlignContent.get_id()
                if sSelfAlignId==sOtherAlignId:
                    tAllHspList.extend(list(oOtherAlignContent.get_hspList()))
                    tAllHspList=sorted(tAllHspList,key=lambda x: x.sbjct_start)
                    bDebug=False
                    #if sSelfAlignId=="ch39_read2381_template_pass_BYK_CB_ONT_1_FAF04998_A":
                        #bDebug=True
                        #print("!!!!!!DEBUG!!!!!!")
                    
                    while True:
                        iBaseHsp=len(tAllHspList)
                        tToRemove=[]
                        tToAdd=[]
                        for iFirstIndex in range(len(tAllHspList)-1):
                            oFirstHsp=tAllHspList[iFirstIndex]
                            iFirstStart=oFirstHsp.query_start
                            iFirstEnd=oFirstHsp.query_end
                            for iSecondIndex in range(iFirstIndex,len(tAllHspList)):
                                oSecondHsp=tAllHspList[iSecondIndex]
                                if oFirstHsp is oSecondHsp:
                                    continue
                                iSecondStart=oSecondHsp.query_start
                                iSecondEnd=oSecondHsp.query_end
                                #print("{},{} : {},{}".format(iFirstStart,iFirstEnd,iSecondStart,iSecondEnd))
                                if iFirstStart>iSecondEnd or iFirstEnd<iSecondStart:
                                    continue
                                #Case 00: Same dual Hsp, keep only one
                                if iFirstStart==iSecondStart and iFirstEnd==iSecondEnd:
                                    #print("Case00")
                                    tToRemove.append(oSecondHsp)
                                    break
                                #Case 01: First Hsp include into Second Hsp, keep only the second
                                if iFirstStart>=iSecondStart and iFirstEnd<=iSecondEnd:
                                    #print("Case01")
                                    tToRemove.append(oFirstHsp)
                                    break
                                #Case 02: First Hsp include the Second Hsp, keep only the first
                                if iFirstStart<=iSecondStart and iFirstEnd>=iSecondEnd:
                                    #print("Case02")
                                    tToRemove.append(oSecondHsp)
                                    break
                                #Case 03: Merge Hsp
                                if (iFirstStart<=iSecondStart and iFirstEnd<=iSecondEnd) or (iFirstStart>=iSecondStart and iFirstEnd>=iSecondEnd):
                                    #print("Case03")
                                    tToRemove.append(oSecondHsp)
                                    tToRemove.append(oFirstHsp)
                                    iQueryStrand=oFirstHsp.frame[0]
                                    sNewQueryPart=None
                                    if iQueryStrand==1:
                                        iNewQueryStart=min(oFirstHsp.query_start,oSecondHsp.query_start)
                                        iNewQueryStop=max(oFirstHsp.query_end,oSecondHsp.query_end)
                                    elif iQueryStrand==-1:
                                        iNewQueryStart=max(oFirstHsp.query_start,oSecondHsp.query_start)
                                        iNewQueryStop=min(oFirstHsp.query_end,oSecondHsp.query_end)
                                    else:
                                        print("ERROR 277 : oHsp strand is not 1 nor -1")
                                    iReadStrand=oFirstHsp.frame[1]
                                    sNewReadPart=None
                                    if iReadStrand==1:
                                        iNewReadStart=min(oFirstHsp.sbjct_start,oSecondHsp.sbjct_start)
                                        iNewReadStop=max(oFirstHsp.sbjct_end,oSecondHsp.sbjct_end)
                                    elif iReadStrand==-1:
                                        iNewReadStart=max(oFirstHsp.sbjct_start,oSecondHsp.sbjct_start)
                                        iNewReadStop=min(oFirstHsp.sbjct_end,oSecondHsp.sbjct_end)
                                    else:
                                        print("ERROR 287 : oHsp strand is not 1 nor -1")
                                    oNewHsp=PseudoBlastHSP(sNewQueryPart,iNewQueryStart,iNewQueryStop,
                                                            sNewReadPart,iNewReadStart,iNewReadStop)
                                    tToAdd.append(oNewHsp)
                                    break
                            if len(tToRemove)!=0:
                                break
                        #if bDebug:print("InitialHspList : {}".format([(X.query_start, X.query_end) for X in tAllHspList]))
                        #if bDebug:print("toRemove : {}".format([(X.query_start, X.query_end) for X in tToRemove]))
                        #if bDebug:print("toAdd : {}".format([(X.query_start, X.query_end) for X in tToAdd]))
                        for oObject in tToRemove:
                            tAllHspList.remove(oObject)
                        for oObject in tToAdd:
                            tAllHspList.extend(tToAdd)
                        tAllHspList=sorted(tAllHspList,key=lambda x: x.sbjct_start)
                        #if bDebug:print("ModifiedHspList : {}".format([(X.query_start, X.query_end) for X in tAllHspList]))
                        
                        #if bDebug:print("-------------Next-------------")
                        if iBaseHsp==len(tAllHspList):
                            break
                    
                    oSelfAlignContent.set_hspList(tAllHspList)
        #Add unexisting AlignContent from the new ToolContent
        setOtherAlignContentId=set([X.get_id() for X in oToolContent.get_alignList()])
        setSelfAlignContentId=set([X.get_id() for X in self.get_alignList()])
        setDeltaOtherAlignContentId=setOtherAlignContentId-setSelfAlignContentId
        for oOtherAlignContent in oToolContent.get_alignList():
            if oOtherAlignContent.get_id() in setDeltaOtherAlignContentId:
                self.add_align(oOtherAlignContent)

class ExonerateContent(ToolContent):
    
    def __init__(self,sOutputFile,sGeneFastaFile,sReadFastaFile,dRead2Strand=None):
        bStrandData=True
        if dRead2Strand is None:
            bStrandData=False
        
        self.strand=None
        self.set_tool("exonerate")
        oBlastMatch=self.exonerate2pseudoBlast(sOutputFile,sGeneFastaFile,sReadFastaFile)
        self.set_queryId(oBlastMatch.query)
        self.set_querySize(oBlastMatch.query_letters)
        self.set_emptyAlignList()
        
        oAlignments=oBlastMatch.alignments
        if len(oAlignments)!=0:
            for oUniqAlignment in oAlignments:
                sRefId=oUniqAlignment.title
                sRefSize=oUniqAlignment.length
                tHspList=oUniqAlignment.hsps
                if len(tHspList)==0:
                    print("WARNING : {} have no Hsp !".format(oUniqAlignment.title))
                    continue
                oAlignContent=AlignContent(sRefId,sRefSize,tHspList)
                if bStrandData:
                    oAlignContent.remove_badStrandHsp(dRead2Strand[sRefId])
                self.add_align(oAlignContent)
                
            oLambdaHsp=tHspList[0]
            iQueryStrand=oLambdaHsp.frame[0]
            self.set_strand(iQueryStrand)
        
    def exonerate2pseudoBlast(self,sFilePath,sGeneFastaFile,sReadFastaFile):
        dDict=self.parse_exonerateFile(sFilePath,sGeneFastaFile,sReadFastaFile)
        dDict=self.refine_exonerateParsing(dDict)
        if len(dDict)!=1:
            exit("ERROR 304 : multiple exonerate ref ??")
        for dbGene in dDict.keys():
            oPseudoBlastMatch=PseudoBlastMatch(dbGene[0],dbGene[1])
            for dbRead in dDict[dbGene]:
                oPseudoBlastAlignment=PseudoBlastAlignments(dbRead[0],dbRead[1])
                for iId in dDict[dbGene][dbRead]:
                    oPseudoBlastHsp=PseudoBlastHSP(
                        #dDict[dbGene][dbRead][iId]["frame"],
                        dDict[dbGene][dbRead][iId]["query"],
                        dDict[dbGene][dbRead][iId]["query_start"],
                        dDict[dbGene][dbRead][iId]["query_end"],
                        dDict[dbGene][dbRead][iId]["sbjct"],
                        dDict[dbGene][dbRead][iId]["sbjct_start"],
                        dDict[dbGene][dbRead][iId]["sbjct_end"]
                        )
                    oPseudoBlastAlignment.hsps.append(oPseudoBlastHsp)
                oPseudoBlastMatch.alignments.append(oPseudoBlastAlignment)
        return oPseudoBlastMatch
    
    def exonerate_reverseComp(self,sString):
        if ">" not in sString and "<" not in sString:
            return "".join([REVERSE_COMP[X] for X in sString[::-1]])
        tString=re.split('[<,>]',sString)
        tNewString=[]
        for iIndex in range(len(tString)-1,-1,-1):
            if "Target Intron" in tString[iIndex]:
                tNewString.append(tString[iIndex])
            else:
                tNewString.append("".join([REVERSE_COMP[X] for X in tString[iIndex][::-1]]))
        sNewString=">".join(tNewString)
        return sNewString
        
        
    def exonerate_reverseAlign(self,sString):
        tString=re.split("\+\+|-\+|\+-|--",sString)
        tResult=[]
        for iIndex in range(len(tString)-1,-1,-1):
            if "|" in tString[iIndex]:
                tResult.append(tString[iIndex][::-1])
            else:
                tResult.append(tString[iIndex])
        sResult="++".join(tResult)
        return sResult
    
    def parse_exonerateFile(self,sFile,sGeneFastaFile,sReadFastaFile):
        dGene2Fasta=Fasta2Dict(sGeneFastaFile)
        if len(dGene2Fasta)!=1:
            exit("ERROR 301 : Many Ref in Ref fasta file")
        dRead2Fasta=Fasta2Dict(sReadFastaFile)
        dResult={}
        bOk=False
        bInAlign=False
        iAlignCounter=0
        sReadAlign=""
        sDataAlign=""
        sGeneAlign=""
        sReadId=""
        sGeneId=""
        iReadStart=-1
        iReadStop=-1
        iGeneStart=-1
        iGeneStop=-1
        iGeneStrand=None
        iReadStrand=None
        iReadSize=None
        iRawScore=None
        for sLine in open(sFile):
            sStripedLine=sLine.strip()
            #print(sStripedLine)
            if "#" in sStripedLine:
                if bInAlign:
                    if iRawScore>=RAWSCORE_THRESHOLD:
                        try:
                            oCrash=dResult[(sGeneId,iGeneSize)]
                        except KeyError:
                            dResult[(sGeneId,iGeneSize)]={}
                        try:
                            oCrash=dResult[(sGeneId,iGeneSize)][(sReadId,iReadSize)]
                        except KeyError:
                            dResult[(sGeneId,iGeneSize)][(sReadId,iReadSize)]={}
                            
                        if iGeneStrand==1:
                            iGeneStart+=1
                        elif iGeneStrand==-1:
                            iGeneStop+=1
                        if iReadStrand==1:
                            iReadStart+=1
                        elif iReadStrand==-1:
                            iReadStop+=1
                                        
                        if iGeneStrand==1:
                            dResult[(sGeneId,iGeneSize)][(sReadId,iReadSize)][len(dResult[(sGeneId,iGeneSize)][(sReadId,iReadSize)])]={
                                        "frame":(iGeneStrand,iReadStrand),
                                        "query":sGeneAlign,
                                        "query_start":iGeneStart-1,
                                        "query_end":iGeneStop-1,
                                        "sbjct":sReadAlign,
                                        "sbjct_start":iReadStart-1,
                                        "sbjct_end":iReadStop-1,
                                        "align_data":sDataAlign,
                                        "raw_score":iRawScore
                                        }
                        else:
                            dResult[(sGeneId,iGeneSize)][(sReadId,iReadSize)][len(dResult[(sGeneId,iGeneSize)][(sReadId,iReadSize)])]={
                                        "frame":(-iGeneStrand,-iReadStrand),
                                        "query":self.exonerate_reverseComp(sGeneAlign),
                                        "query_start":iGeneStop-1,
                                        "query_end":iGeneStart-1,
                                        "sbjct":self.exonerate_reverseComp(sReadAlign),
                                        "sbjct_start":iReadStop-1,
                                        "sbjct_end":iReadStart-1,
                                        "align_data":self.exonerate_reverseAlign(sDataAlign),
                                        "raw_score":iRawScore
                                        }
                    #print("sReadId",sReadId,iReadStart,iReadStop)
                    #print("sGeneId",sGeneId,iGeneStart,iGeneStop)
                    #print(sReadAlign)
                    #print(sDataAlign)
                    #print(sGeneAlign)
                    #exit("DONE")                
                bOk=False
                bInAlign=False
                iAlignCounter=0
                sReadAlign=""
                sDataAlign=""
                sGeneAlign=""
                sReadId=""
                sGeneId=""
                iReadStart=-1
                iReadStop=-1
                iGeneStart=-1
                iGeneStop=-1
                iGeneStrand=None
                iReadStrand=None
                iReadSize=None
                iRawScore=None
                continue
            if "C4 Alignment" in sStripedLine:
                bOk=True
                continue
            if bOk:
                if "Query:" in sStripedLine:
                    tLine=sStripedLine.split(": ")
                    sReadId=tLine[-1]
                    try:
                        iReadSize=len(dRead2Fasta[sReadId])
                    except KeyError:
                        exit("ERROR 433 : Read {} in exonerate file is not present in Read fasta".format(sReadId))
                    iReadStrand=1
                    #print(sReadId)
                elif "Target:" in sStripedLine:
                    if " [revcomp]" in sStripedLine:
                        sStripedLine=sStripedLine.replace(" [revcomp]","")
                        iGeneStrand=-1
                    else:
                        iGeneStrand=1
                    tLine=sStripedLine.split(": ")                    
                    sGeneId=tLine[-1]
                    try:
                        iGeneSize=len(dGene2Fasta[sGeneId])
                    except KeyError:
                        exit("ERROR 447 : Ref in fasta file do not correpond to Ref in exonerate file")
                    #print(sGeneId)
                elif "Raw score:" in sStripedLine:
                    tLine=sStripedLine.split(": ")
                    iRawScore=int(tLine[-1])
                elif "Query range:" in sStripedLine:
                    tLine=sStripedLine.split(": ")
                    tItems=tLine[-1].split("->")
                    iReadStart=int(tItems[0])
                    iReadStop=int(tItems[-1])
                    #print(iReadStart,iReadStop)
                elif "Target range:" in sStripedLine:
                    tLine=sStripedLine.split(": ")
                    tItems=tLine[-1].split("->")
                    iGeneStart=int(tItems[0])
                    iGeneStop=int(tItems[-1])
                    bInAlign=True
                    iAlignCounter=0
                    #print(iGeneStart,iGeneStop)
                elif bInAlign:
                    if len(sStripedLine)==0:
                        continue
                    if iAlignCounter==0:
                        iFirstLim=sLine.find(": ")+2
                        iLastLim=iFirstLim+sLine[iFirstLim:].find(" :")
                        sSeq=sLine[iFirstLim:iLastLim]
                        sReadAlign+=sSeq
                        iAlignCounter+=1
                    elif iAlignCounter==2:
                        sSeq=sLine[iFirstLim:iLastLim]
                        sGeneAlign+=sSeq
                        iAlignCounter=0
                    elif iAlignCounter==1:
                        sSeq=sLine[iFirstLim:iLastLim]
                        sDataAlign+=sSeq
                        iAlignCounter+=1
        return dResult
        
    def refine_exonerateParsing(self,dDict):
        for dbGeneKey in dDict.keys():
            for dbReadKey in sorted(dDict[dbGeneKey]):
                #print(dbReadKey)
                iAlignId=-1
                #print(dDict[dbGeneKey][dbReadKey])
                while iAlignId<max(dDict[dbGeneKey][dbReadKey]):
                    iAlignId+=1
                    iChevronSpace=2
                    iChevronSize=4
                    tChevron=["<<<<",">>>>"]
                    #print(iAlignId)
                    #print(dDict[dbGeneKey][dbReadKey].keys())
                    while "Intron" in dDict[dbGeneKey][dbReadKey][iAlignId]["sbjct"]:
                        sCurrentDataAlign=dDict[dbGeneKey][dbReadKey][iAlignId]["align_data"]
                        sNewDataAlign=""
                        sCurrentGeneSeq=dDict[dbGeneKey][dbReadKey][iAlignId]["query"]
                        sNewGeneSeq=""
                        iCurrentGeneStart=dDict[dbGeneKey][dbReadKey][iAlignId]["query_start"]
                        iNewGeneStart=-1
                        iCurrentGeneStop=dDict[dbGeneKey][dbReadKey][iAlignId]["query_end"]
                        iNewGeneStop=-1
                        sCurrentReadSeq=dDict[dbGeneKey][dbReadKey][iAlignId]["sbjct"]
                        sNewReadSeq=""
                        iCurrentReadStart=dDict[dbGeneKey][dbReadKey][iAlignId]["sbjct_start"]
                        iNewReadStart=-1
                        iCurrentReadStop=dDict[dbGeneKey][dbReadKey][iAlignId]["sbjct_end"]
                        iNewGeneStop=-1
                        
                        ##Find Intron Lim in the Alignment
                        iLimInf=max(sCurrentReadSeq.find(tChevron[0])-iChevronSpace,sCurrentReadSeq.find(tChevron[1])-iChevronSpace)
                        iLimSup=max(iLimInf+iChevronSize*2+sCurrentReadSeq[iLimInf+iChevronSize:].find(tChevron[0])+iChevronSpace,
                                iLimInf+iChevronSize*2+sCurrentReadSeq[iLimInf+iChevronSize:].find(tChevron[1])+iChevronSpace)
                        iIntronSize=[int(s) for s in sCurrentDataAlign[iLimInf:iLimSup].split() if s.isdigit()][0]
                        
                        ##Update data
                        sNewGeneSeq=sCurrentGeneSeq[iLimSup:]
                        sCurrentGeneSeq=sCurrentGeneSeq[:iLimInf]
                        iCurrentGeneStart=iCurrentGeneStart
                        iNewGeneStop=iCurrentGeneStop
                        if dDict[dbGeneKey][dbReadKey][iAlignId]["frame"][0]==1:
                            iCurrentGeneStop=iCurrentGeneStart+len(sCurrentGeneSeq.replace("-",""))-1
                            iNewGeneStart=iCurrentGeneStop+iIntronSize+1
                        elif dDict[dbGeneKey][dbReadKey][iAlignId]["frame"][0]==-1:
                            iCurrentGeneStop=iCurrentGeneStart-len(sCurrentGeneSeq.replace("-",""))+1
                            iNewGeneStart=iCurrentGeneStop-iIntronSize-1
                        
                        sNewReadSeq=sCurrentReadSeq[iLimSup:]
                        sCurrentReadSeq=sCurrentReadSeq[:iLimInf]
                        iCurrentReadStart=iCurrentReadStart
                        iNewReadStop=iCurrentReadStop
                        if dDict[dbGeneKey][dbReadKey][iAlignId]["frame"][-1]==1:
                            iCurrentReadStop=iCurrentReadStart+len(sCurrentReadSeq.replace("-",""))-1
                            iNewReadStart=iCurrentReadStop+1
                        elif dDict[dbGeneKey][dbReadKey][iAlignId]["frame"][-1]==-1:
                            iCurrentReadStop=iCurrentReadStart-len(sCurrentReadSeq.replace("-",""))+1
                            iNewReadStart=iCurrentReadStop-1
                        
                        sNewDataAlign=sCurrentDataAlign[iLimSup:]
                        sCurrentDataAlign=sCurrentDataAlign[:iLimInf]
                        
                        ##Update current data
                        dDict[dbGeneKey][dbReadKey][iAlignId]={
                                    "frame":dDict[dbGeneKey][dbReadKey][iAlignId]["frame"],
                                    "query":sCurrentGeneSeq,
                                    "query_start":iCurrentGeneStart,
                                    "query_end":iCurrentGeneStop,
                                    "sbjct":sCurrentReadSeq,
                                    "sbjct_start":iCurrentReadStart,
                                    "sbjct_end":iCurrentReadStop,
                                    "align_data":sCurrentDataAlign
                                    }
                        
                        ##Add new data
                        iDebug=len(dDict[dbGeneKey][dbReadKey])
                        dDict[dbGeneKey][dbReadKey][len(dDict[dbGeneKey][dbReadKey])]={
                                    "frame":dDict[dbGeneKey][dbReadKey][iAlignId]["frame"],
                                    "query":sNewGeneSeq,
                                    "query_start":iNewGeneStart,
                                    "query_end":iNewGeneStop,
                                    "sbjct":sNewReadSeq,
                                    "sbjct_start":iNewReadStart,
                                    "sbjct_end":iNewReadStop,
                                    "align_data":sNewDataAlign
                                    }
                                    
                        if DEBUG_BOOL and dbReadKey[0]==DEBUG_READID:
                            print("Current:",dDict[dbGeneKey][dbReadKey][iAlignId]["frame"])
                            print("Next:",dDict[dbGeneKey][dbReadKey][iDebug]["frame"])
        return dDict

class YassContent(ToolContent):
    
    def __init__(self,sFile):
        self.strand=None
        self.set_tool("yass")
        oBlastMatch=self.yass2pseudoBlast(sFile)
        self.set_queryId(oBlastMatch.query)
        self.set_querySize(oBlastMatch.query_letters)
        self.set_emptyAlignList()
        
        oAlignments=oBlastMatch.alignments
        if len(oAlignments)!=0:
            for oUniqAlignment in oAlignments:
                sRefId=oUniqAlignment.title
                sRefSize=oUniqAlignment.length
                tHspList=oUniqAlignment.hsps
                oAlignContent=AlignContent(sRefId,sRefSize,tHspList)
                self.add_align(oAlignContent)
                
            oLambdaHsp=tHspList[0]
            iQueryStrand=oLambdaHsp.frame[0]
            self.set_strand(iQueryStrand)
    
    def yass2pseudoBlast(self,sFilePath):
        dDict=self.parse_yassFile(sFilePath)
        if len(dDict)!=1:
            exit("ERROR 141 : multiple yass ref ??")
        for dbGene in dDict.keys():
            oPseudoBlastMatch=PseudoBlastMatch(dbGene[0],dbGene[1])
            for dbRead in dDict[dbGene]:
                oPseudoBlastAlignment=PseudoBlastAlignments(dbRead[0],dbRead[1])
                for iId in dDict[dbGene][dbRead]:
                    oPseudoBlastHsp=PseudoBlastHSP(
                        #dDict[dbGene][dbRead][iId]["frame"],
                        dDict[dbGene][dbRead][iId]["query"],
                        dDict[dbGene][dbRead][iId]["query_start"],
                        dDict[dbGene][dbRead][iId]["query_end"],
                        dDict[dbGene][dbRead][iId]["sbjct"],
                        dDict[dbGene][dbRead][iId]["sbjct_start"],
                        dDict[dbGene][dbRead][iId]["sbjct_end"]
                        )
                    oPseudoBlastAlignment.hsps.append(oPseudoBlastHsp)
                oPseudoBlastMatch.alignments.append(oPseudoBlastAlignment)
        return oPseudoBlastMatch
        
    def parse_yassFile(self,sFilePath):    
        dContent={}
        iAsterisk=0
        sGenePartContent=""
        sReadPartContent=""
        for sLine in open(sFilePath):
            sLine=sLine.strip()
            if len(sLine)==0:
                iBlank=0
                continue
            elif sLine[0]=="*" :
                iAsterisk+=1
                if iAsterisk==5:
                    dContent[dbGene][dbRead][len(dContent[dbGene][dbRead])]={
                        "frame":None,
                        "query":sGenePartContent,
                        "query_start":iGenePartStart,
                        "query_end":iGenePartEnd,
                        "sbjct":sReadPartContent,
                        "sbjct_start":iReadPartStart,
                        "sbjct_end":iReadPartEnd
                        }
                    iAsterisk=1
                    sGenePartContent=""
                    sReadPartContent=""
                if iAsterisk==1:
                    sDelimiter=re.compile(r"[-\(\)]")
                    tLine=sDelimiter.split(sLine)
                    
                    if int(tLine[1])<int(tLine[2]):
                        iGenePartStart=int(tLine[1])
                        iGenePartEnd=int(tLine[2])
                        iReadPartStart=int(tLine[4])
                        iReadPartEnd=int(tLine[5])
                    else:
                        iGenePartStart=int(tLine[2])
                        iGenePartEnd=int(tLine[1])
                        iReadPartStart=int(tLine[5])
                        iReadPartEnd=int(tLine[4])
                if iAsterisk==2:
                    tLine=sLine.split("/")
                    sGeneLine=tLine[0]
                    sReadLine=tLine[1]
                    sGeneId=sGeneLine.split('"')[1]
                    iGeneSize=int(sGeneLine.split("(")[-1].split(" ")[0])
                    sReadId=sReadLine.split('"')[1]
                    iReadSize=int(sReadLine.split("(")[-1].split(" ")[0])
                    
                    dbGene=tuple([sGeneId,iGeneSize])
                    dbRead=tuple([sReadId,iReadSize])
                    try:
                        oCrash=dContent[dbGene]
                        try:
                            oCrash=dContent[dbGene][dbRead]
                        except KeyError:
                            dContent[dbGene][dbRead]={}
                    except KeyError:
                        dContent[dbGene]={}
                        dContent[dbGene][dbRead]={}
            else:
                iBlank+=1
                if iBlank==2:
                    sGenePartContent+=sLine
                if iBlank==4:
                    sReadPartContent+=sLine
        dContent[dbGene][dbRead][len(dContent[dbGene][dbRead])]={
            "frame":None,
            "query":sGenePartContent,
            "query_start":iGenePartStart,
            "query_end":iGenePartEnd,
            "sbjct":sReadPartContent,
            "sbjct_start":iReadPartStart,
            "sbjct_end":iReadPartEnd
            }
        return dContent

class BlastContent(ToolContent): 
       
    def __init__(self,sFile):
        self.strand=None
        self.set_tool("blast")
        oFileContent=open(sFile)
        oAllBlastMatchs=NCBIXML.parse(oFileContent)
        bFirst=True
        for oBlastMatch in oAllBlastMatchs:
            oAlignments=oBlastMatch.alignments
            if bFirst:
                bFirst=False
            else:
                exit("ERROR 48 : multiple blast match ??")
            self.set_queryId(oBlastMatch.query)
            self.set_querySize(oBlastMatch.query_letters)
            self.set_emptyAlignList()
            
            if len(oAlignments)!=0:
                for oUniqAlignment in oAlignments:
                    sRefId=oUniqAlignment.title
                    sRefId=sRefId.replace(" No definition line","")
                    sRefSize=oUniqAlignment.length
                    tHspList=oUniqAlignment.hsps
                    oAlignContent=AlignContent(sRefId,sRefSize,tHspList)
                    self.add_align(oAlignContent)
                    
                oLambdaHsp=tHspList[0]
                iQueryStrand=oLambdaHsp.frame[0]
                self.set_strand(iQueryStrand)

class AlignContent:
    
    def __init__(self,sId,sSize,tObjectList):
        self.align_id=sId
        self.align_size=sSize
        self.align_hsp=tObjectList
        self.exonCovering={}
        
        self.cumul_read_strand=None
        self.cumul_oriented_read_size=None
        self.cumul_gene_strand=None
        self.cumul_oriented_gene_size=None
        self.assign_cumul_strand()
        
        self.read_strand=None
        self.gene_strand=None
        self.assign_major_strand()

    def get_read_size(self):
        return self.align_size
    
    def remove_badStrandHsp(self,dStrand):
        tToRemove=[]
        iMegablastGeneStrand=dStrand["geneStrand"]
        iMegablastReadStrand=dStrand["readStrand"]
        for oHsp in self.get_hspList():
            if oHsp.frame[0]==iMegablastGeneStrand:
                if oHsp.frame[1]!=iMegablastReadStrand:
                    tToRemove.append(oHsp)
            elif oHsp.frame[0]==-iMegablastGeneStrand:
                if oHsp.frame[1]==iMegablastReadStrand:
                    tToRemove.append(oHsp)
        tNewHspList=[X for X in self.get_hspList() if X not in tToRemove]
        self.set_hspList(tNewHspList)
    
    def assign_cumul_strand(self):
        iGeneStrand=0
        iReadStrand=0
        #iOrientedGeneSize=0
        #iOrientedReadSize=0
        #print("========================")
        for oHsp in self.get_hspList():
            #print(oHsp.frame)
            iGeneStrand+=oHsp.frame[0]
            iReadStrand+=oHsp.frame[1]
            #iOrientedGeneSize+=len(oHsp.query)*oHsp.frame[0]
            #iOrientedReadSize+=len(oHsp.sbjct)*oHsp.frame[1]
        #print("---------------")
        #print(iGeneStrand,iReadStrand)
        #print(iOrientedGeneSize,iOrientedReadSize)
        self.set_cumul_gene_strand(iGeneStrand)
        self.set_cumul_read_strand(iReadStrand)
        #self.set_cumul_oriented_read_size(iOrientedGeneSize)
        #self.set_cumul_oriented_gene_size(iOrientedReadSize)
    
    #def set_cumul_oriented_read_size(self,iInt):
        #self.cumul_oriented_read_size=iInt
    
    #def set_cumul_oriented_gene_size(self,iInt):
        #self.cumul_oriented_gene_size=iInt
    
    #def get_cumul_oriented_read_size(self):
        #return self.cumul_oriented_read_size
    
    #def get_cumul_oriented_gene_size(self):
        #return self.cumul_oriented_gene_size
    
    def assign_major_strand(self):
        if self.get_cumul_gene_strand()>0:
            self.set_gene_strand(1)
        elif self.get_cumul_gene_strand()<0:
            self.set_gene_strand(-1)
        if self.get_cumul_read_strand()>0:
            self.set_read_strand(1)
        elif self.get_cumul_read_strand()<0:
            self.set_read_strand(-1)
    
    def set_gene_strand(self,iInt):
        self.gene_strand=iInt
    
    def get_gene_strand(self):
        return self.gene_strand
    
    def set_read_strand(self,iInt):
        self.read_strand=iInt
    
    def get_read_strand(self):
        return self.read_strand
    
    def set_cumul_read_strand(self,iInt):
        self.cumul_read_strand=iInt
    
    def get_cumul_read_strand(self):
        return self.cumul_read_strand
    
    def set_cumul_gene_strand(self,iInt):
        self.cumul_gene_strand=iInt
    
    def get_cumul_gene_strand(self):
        return self.cumul_gene_strand
    
    def get_id(self):
        return self.align_id
    
    def get_hspList(self):
        return self.align_hsp
    
    def set_hspList(self,tObjectList):
        self.align_hsp=list(tObjectList)
    
    def get_alignContent_size(self):
        return self.align_size
        
    def set_alignContent_selfCover(self,fFloat):
        self.align_self_cover=fFloat
        
    def get_alignContent_selfCover(self):
        return self.align_self_cover
    
    def compute_ref_cover(self,iInt):
        iDelta=max(self.get_globalEnd(),self.get_globalStart())-min(self.get_globalEnd(),self.get_globalStart())
        fRatio=float(iDelta)/iInt*100
        self.set_alignContent_refCover(fRatio)

    def set_alignContent_refCover(self,fFloat):
        self.align_ref_cover=fFloat
        
    def get_alignContent_refCover(self):
        return self.align_ref_cover
        
    def compute_selfGap(self):
        dbValue=self.get_globalSelfPosition()
        iGap=0
        iCover=0
        for iIndex in range(len(dbValue)-1):
            iGap+=1
            dbOne=dbValue[iIndex]
            dbTwo=dbValue[iIndex+1]
            iStop=dbOne[-1]+1
            iStart=dbTwo[0]
            iDelta=iStart-iStop
            iCover+=iDelta
        fCover= float(iCover)/self.get_alignContent_size()*100
        self.set_number_of_selfGap(iGap)
        self.set_selfGap_covering(fCover)
        
    def set_number_of_selfGap(self,iInt):
        self.nb_self_gap=iInt
        
    def get_number_of_selfGap(self):
        return self.nb_self_gap
    
    def set_selfGap_covering(self,fFloat):
        self.self_gap_covering=fFloat
        
    def get_selfGap_covering(self):
        return self.self_gap_covering
        
    def get_selfStrand(self):
        return self.selfStrand
        
    def set_selfStrand(self,iInt):
        self.selfStrand=iInt
        
    def get_refStrand(self):
        return self.refStrand
        
    def set_refStrand(self,iInt):
        self.refStrand=iInt
    
    def compute_self_cover(self):
        iSelfSize=self.get_alignContent_size()
        bPrevious=False
        iPreviousSelfStop=None
        iSumSelfOverlap=0
        iSelfCover=0
        
        #print(self.align_id)
        
        #DEBUG
        iStrand=None
        iRefStrand=None
        for oHsp in self.get_hspList():
            if iStrand==None:
                iStrand=oHsp.frame[1] #Read strand
                iRefStrand=oHsp.frame[0]
            elif iStrand!=oHsp.frame[1] or iRefStrand!=oHsp.frame[0]:
                print("WARNING 499 : divergente strand for the same read alignment")
                print(self.get_id(),self.get_alignContent_size())
                for oHsp in self.get_hspList():
                    print(oHsp.frame,oHsp.query_start,oHsp.query_end,oHsp.sbjct_start,oHsp.sbjct_end)
                print("selfCover arbitrary set at -1")
                self.set_alignContent_selfCover(-1.0)
                return None
        
        self.set_refStrand(oHsp.frame[0])
        self.set_selfStrand(oHsp.frame[1])
        #print("--check--",self.get_refStrand(),self.get_selfStrand(),self.get_alignContent_size())
        
        #print(iSelfSize,iStrand)
        
        if iStrand==1:        
            for oHsp in sorted(self.get_hspList(),key=lambda x: x.sbjct_start):
                iSelfStart=oHsp.sbjct_start
                iSelfStop=oHsp.sbjct_end
                iSelfRealSize=iSelfStop-iSelfStart+1
                #print(iSelfStart,iSelfStop,iSelfRealSize,oHsp.frame,oHsp.query_start,oHsp.query_)                
                if bPrevious:
                    if iSelfStart<=iPreviousSelfStop+1:
                        iSumSelfOverlap=iPreviousSelfStop-iSelfStart+1
                    
                bPrevious=True
                iPreviousSelfStop=iSelfStop
                iSelfCover+=iSelfRealSize
                #print(iSelfCover,iSumSelfOverlap)
        elif iStrand==-1:
            for oHsp in sorted(self.get_hspList(),key=lambda x: x.sbjct_end):
                iSelfStart=oHsp.sbjct_end
                iSelfStop=oHsp.sbjct_start
                iSelfRealSize=iSelfStop-iSelfStart+1
                #print(iSelfStart,iSelfStop,iSelfRealSize,oHsp.frame)                
                if bPrevious:
                    if iSelfStart<=iPreviousSelfStop+1:
                        iSumSelfOverlap=iPreviousSelfStop-iSelfStart+1
                    
                bPrevious=True
                iPreviousSelfStop=iSelfStop
                iSelfCover+=iSelfRealSize
        else:
            exit("ERROR 523 : strand in Hsp is not 1 nor -1")
                
        iSelfCover=iSelfCover-iSumSelfOverlap
        fSelfCover=float(iSelfCover)/iSelfSize*100
        self.set_alignContent_selfCover(fSelfCover)
        
    def define_position(self):
        """
        Scan all Query/Self position in order to fusion element with overlap
        Pre_requisites for assign self_exonCovering
        """
        tAllQueryPosition=[] #Gene/ref
        tAllSelfPosition=[] #sbjct/reads
        
        if self.get_alignContent_selfCover()==-1.0:
            ##Bad items. Assign arbitrary value
            self.set_globalStart(-1)
            self.set_globalEnd(-1)
            self.set_number_of_selfGap(-1)
            self.set_selfGap_covering(-1)
            return None
        
        #print(self.get_alignContent_selfCover())
        #print(self.get_alignContent_size())
        
        if self.get_refStrand()==1:
            for oHsp in sorted(self.get_hspList(),key=lambda x: x.query_start):
                iQueryStart=oHsp.query_start
                iQueryStop=oHsp.query_end
                tAllQueryPosition.append((iQueryStart,iQueryStop))
        else:
            for oHsp in sorted(self.get_hspList(),key=lambda x: x.query_end):
                iQueryStart=oHsp.query_end
                iQueryStop=oHsp.query_start
                tAllQueryPosition.append((iQueryStart,iQueryStop))
        
        if self.get_selfStrand()==1:
            #print([(X.sbjct_start,X.sbjct_end) for X in sorted(self.get_hspList(),key=lambda x: x.sbjct_start)])
            for oHsp in sorted(self.get_hspList(),key=lambda x: x.sbjct_start):    
                iSelfStart=oHsp.sbjct_start
                iSelfStop=oHsp.sbjct_end
                tAllSelfPosition.append((iSelfStart,iSelfStop))
        else:
            #print([(X.sbjct_start,X.sbjct_end) for X in sorted(self.get_hspList(),key=lambda x: x.sbjct_end)])
            for oHsp in sorted(self.get_hspList(),key=lambda x: x.sbjct_end):    
                iSelfStart=oHsp.sbjct_end
                iSelfStop=oHsp.sbjct_start
                tAllSelfPosition.append((iSelfStart,iSelfStop))
        #print(tAllSelfPosition)
        
        tGlobalQueryPosition=[]
        for dbQueryPosition in tAllQueryPosition:
            if len(tGlobalQueryPosition)==0:
                tGlobalQueryPosition.append(dbQueryPosition)
            elif dbQueryPosition[0]<tGlobalQueryPosition[-1][1]:
                #Fusion
                tGlobalQueryPosition[-1]=(min(dbQueryPosition[0],tGlobalQueryPosition[-1][0]),max(dbQueryPosition[1],tGlobalQueryPosition[-1][1]))
            else:
                tGlobalQueryPosition.append(dbQueryPosition)
                
        tGlobalSelfPosition=[]
        for dbSelfPosition in tAllSelfPosition:
            if len(tGlobalSelfPosition)==0:
                tGlobalSelfPosition.append(dbSelfPosition)
            elif dbSelfPosition[0]<tGlobalSelfPosition[-1][1]:
                #Fusion
                tGlobalSelfPosition[-1]=(min(dbSelfPosition[0],tGlobalSelfPosition[-1][0]),max(dbSelfPosition[1],tGlobalSelfPosition[-1][1]))
            else:
                tGlobalSelfPosition.append(dbSelfPosition)
        
        dbGlobalQueryPosition=tuple(tGlobalQueryPosition)
        self.set_globalQueryPosition(dbGlobalQueryPosition)
        dbGlobalSelfPosition=tuple(tGlobalSelfPosition)
        self.set_globalSelfPosition(dbGlobalSelfPosition)
        
        self.set_globalStart(dbGlobalQueryPosition[0][0])
        self.set_globalEnd(dbGlobalQueryPosition[-1][-1])
        
        self.compute_selfGap()
        
    def set_globalStart(self,iValue):
        self.globalStart=iValue
        
    def get_globalStart(self):
        return self.globalStart
    
    def set_globalEnd(self,iValue):
        self.globalEnd=iValue
        
    def get_globalEnd(self):
        return self.globalEnd
            
    def set_globalQueryPosition(self,dbTuple):
        self.globalQueryPosition=dbTuple
    
    def get_globalQueryPosition(self):
        return self.globalQueryPosition
        
    def set_globalSelfPosition(self,dbTuple):
        self.globalSelfPosition=dbTuple
        
    def get_globalSelfPosition(self):
        return self.globalSelfPosition
        
    def describe_self(self,oTrContent=None,bStdout=True):        
        sContent=""
        dbTarget=(self.get_gene_strand(),self.get_read_strand(),
                    self.get_id(),self.get_alignContent_size(),
                    self.get_alignContent_selfCover(),
                    #self.get_globalStart(),self.get_globalEnd()
                    self.get_selfGap_covering(),self.get_number_of_selfGap(),
                    self.get_alignContent_refCover())
        for oValue in dbTarget:
            if sContent!="":
                sContent+="\t"
            sContent+="{}".format(oValue)
        dTarget=self.get_exonCovering()
        iSumPart=0
        if oTrContent is None:
            for sKey in sorted(dTarget):
                if dTarget[sKey]>0:
                    iSumPart+=1
                sContent+="\t{}".format(dTarget[sKey]["deltaStart"])
                sContent+="\t{}".format(dTarget[sKey]["ExonCovering"])
                sContent+="\t{}".format(dTarget[sKey]["deltaStop"])
        else:
            #print("-------",self.get_alignContent_size(),"-------")
            for oExonContent in oTrContent.get_exonList():
                sKey=oExonContent.get_index()
                if dTarget[sKey]["ExonCovering"]>0:
                    iSumPart+=1
                sContent+="\t{}".format(dTarget[sKey]["deltaStart"])
                sContent+="\t{}".format(dTarget[sKey]["ExonCovering"])
                sContent+="\t{}".format(dTarget[sKey]["deltaStop"])
                
        sContent+="\t{}".format(iSumPart)
        if bStdout:
            print(sContent)
        sContent+="\n"
        return sContent
        
    def describe_header(self,oTrContent=None,bStdout=True):
        #sContent="Tool\tReadId\tReadSize\t%ReadCover\t%ReadGap\t#Gap\t%GapedGeneCover"
        sContent="Tool\tGeneStrand\tReadStrand\tReadId\tReadSize\t%ReadCover\t%ReadGap\t#Gap\t%GapedGeneCover"
        if oTrContent is None:
            for sKey in sorted(self.get_exonCovering()):
                sContent+="\t{}.e5".format(sKey)
                sContent+="\t{}".format(sKey)
                sContent+="\t{}.e3".format(sKey)
        else:
            for oExonContent in oTrContent.get_exonList():
                sContent+="\t{}.e5".format(oExonContent.get_index().split(".")[0])
                sContent+="\t{}".format(oExonContent.get_index())
                sContent+="\t{}.e3".format(oExonContent.get_index().split(".")[0])
        sContent+="\tSumExon"
        if bStdout:
            print(sContent)
        sContent+="\n"
        return sContent
    
    def describe_empty(self,bStdout=True):
        sContent=""
        dbTarget=(self.get_gene_strand(),self.get_read_strand(),
                    self.get_id(),self.get_alignContent_size(),
                    self.get_alignContent_selfCover(),
                    #self.get_globalStart(),self.get_globalEnd()
                    self.get_selfGap_covering(),self.get_number_of_selfGap(),
                    self.get_alignContent_refCover())
        sContent+="{}".format(self.get_id())
        for oValue in dbTarget[1:]:
            if sContent!="":
                sContent+="\t"
            sContent+=""
        sContent+="\t"
        if bStdout:
            print(sContent)
        sContent+="\n"
        return sContent
    
    def assign_exonCovering(self,oTranscriptContent):
        tExonList=oTranscriptContent.get_exonList()
        dExId2Include={}

        for oExon in tExonList:
            #dExId2Include[oExon.get_index()]=0.00
            dExId2Include[oExon.get_index()]={"ExonCovering":0.00,
                                                "deltaStart":0,
                                                "deltaStop":0}
            
        if self.get_alignContent_selfCover()==-1.0:
            self.set_exonCovering(dExId2Include)
            return None
        
        for dbPosition in self.get_globalQueryPosition():
            iCurrentStart=dbPosition[0]
            iCurrentStop=dbPosition[1]
                        
            for oExon in tExonList:
                sExId=oExon.get_index()
                iExStart=oExon.get_start()
                iExStop=oExon.get_stop()
                iExLength=iExStop-iExStart+1
                
                #print(iCurrentStart,iCurrentStop,"vs",iExStart,iExStop)
                #print("'-> DeltaValue",iExStart,iExStop)
                #print("'-> StoredValue",dExId2Include[sExId]["deltaStart"],dExId2Include[sExId]["deltaStop"])
                if iCurrentStart>iExStop or iCurrentStop<iExStart:
                    continue
                else: #Hsp cover/is include/overlap the exon
                    iDeltaStart=iExStart-iCurrentStart
                    iDeltaStop=iCurrentStop-iExStop
                    fCoveringValue=0.0
                    #Case 01: FusionnedHsp is include into the exon
                    if iCurrentStart>=iExStart and iCurrentStop<=iExStop:
                        fCoveringValue=round((iCurrentStop-iCurrentStart+1)/(iExLength)*100,2)
                    #Case 02: FusionnedHsp cover all the exon
                    if iCurrentStart<=iExStart and iCurrentStop>=iExStop:
                        fCoveringValue=100.00
                    #Case 03: FusionnedHsp overlap the start of the exon
                    if iCurrentStart<=iExStart and iCurrentStop<=iExStop:
                        fCoveringValue=round((iCurrentStop-iExStart+1)/(iExLength)*100,2)
                    #Case 04: FusionnedHsp overlap the end of the exon
                    if iCurrentStart>=iExStart and iCurrentStop>=iExStop:
                        fCoveringValue=round((iExStop-iCurrentStart+1)/(iExLength)*100,2)
                    #dExId2Include[sExId]+=fCoveringValue
                    dExId2Include[sExId]["ExonCovering"]+=fCoveringValue
                    #print(dExId2Include[sExId]["deltaStart"],dExId2Include[sExId]["deltaStop"],
                            #dExId2Include[sExId]["deltaStart"]!=0,dExId2Include[sExId]["deltaStop"]!=0,
                            #dExId2Include[sExId]["deltaStart"]!=0 or dExId2Include[sExId]["deltaStop"]!=0)
                    if dExId2Include[sExId]["deltaStart"]==0 or dExId2Include[sExId]["deltaStop"]==0:
                        dExId2Include[sExId]["deltaStart"]=iDeltaStart
                        dExId2Include[sExId]["deltaStop"]=iDeltaStop
                        #print("delta",iDeltaStart,iDeltaStop)
                        #print("-------------")
                    else:
                        """Data read from left to right. If different Hsp for an exon, update onlyt the deltaStop value"""
                        #exit("ERROR 786 : multiple non fusionned Hsp cover the same exon")
                        dExId2Include[sExId]["deltaStop"]=iDeltaStop
        #print(dExId2Include)

        self.set_exonCovering(dExId2Include)
            
    def set_exonCovering(self,dDict):
        self.exon_covering=dDict
    
    def get_exonCovering(self):
        return self.exon_covering

class TranscriptContent:
    
    def __init__(self,sGeneId,sRefFile):
        self.gene_start=None
        self.gene_end=None
        self.gene_strand=None
        self.transcript_to_exon={}
        self.set_emptyExonList()
        
        self.parse_gffFile(sRefFile,sGeneId)
        self.fusion_firstCommonStart()
        self.fusion_lastCommonStop()
        self.sort_exonList()
        self.assign_exonId()
        self.update_transcriptList()
        self.correct_exonCoord()
    
    def parse_gffFile(self,sPathFile,sGeneId):
        sTrId=None
        
        for sLine in open(sPathFile):
            sLine=sLine.replace("\n","")
            tLine=sLine.split("\t")
            #tRef=["seqname","source","feature","start","end","score","strand","frame","attribute"]
            dRef2Content=dict(zip(GFF_COLUMN,tLine))
            tAttribute=dRef2Content["attribute"].split(";")
            sCurrentGeneId=tAttribute[-1].split(" ")[-1]
            if len(tAttribute)>1:
                sCurrentTrId=tAttribute[-2].split(" ")[-1]
            if len(tAttribute)>2:
                sCurrentExonId=tAttribute[-3].split(" ")[-1]
            if sCurrentGeneId==sGeneId and dRef2Content["feature"].lower()==FEATURE_TRANSCRIPT:
                sTrId=sCurrentTrId
                self.add_transcript(sTrId)
            elif sCurrentGeneId==sGeneId and dRef2Content["feature"].lower()==FEATURE_EXON and sTrId in dRef2Content["attribute"]:
                oExon=ExonContent(int(dRef2Content["start"]),int(dRef2Content["end"]),sTrId)
                self.update_exonList(oExon)
            elif sCurrentGeneId==sGeneId and dRef2Content["feature"].lower()==FEATURE_GENE:
                self.set_geneCoord(int(dRef2Content["start"]),int(dRef2Content["end"]),STRAND_TO_INT[dRef2Content["strand"]])
    
    def fusion_lastCommonStop(self):
        self.sort_exonList()
        iLastStart=None
        tExonCommonStop=[]
        for oExon in self.get_exonList()[::-1]:
            if not iLastStart:
                iLastStart=oExon.get_start()
                tExonCommonStop.append(oExon)
            elif iLastStart==oExon.get_start():
                tExonCommonStop.append(oExon)
        if len(tExonCommonStop)>1:
            oRefExon=tExonCommonStop[0]
            for oExon in tExonCommonStop[1:]:
                oRefExon.update(oExon)
                self.remove_exon(oExon)
    
    def fusion_firstCommonStart(self):
        self.sort_exonList()
        iFirstStop=None
        tExonCommonStart=[]
        for oExon in self.get_exonList():
            if not iFirstStop:
                iFirstStop=oExon.get_stop()
                tExonCommonStart.append(oExon)
            elif iFirstStop==oExon.get_stop():
                tExonCommonStart.append(oExon)
        if len(tExonCommonStart)>1:
            oRefExon=tExonCommonStart[0]
            for oExon in tExonCommonStart[1:]:
                oRefExon.update(oExon)
                self.remove_exon(oExon)
    
    def set_geneCoord(self,iStart,iEnd,iStrand):
        self.gene_start=iStart
        self.gene_end=iEnd
        self.gene_strand=iStrand
        
    def get_geneStart(self):
        return self.gene_start
        
    def get_geneEnd(self):
        return self.gene_end
        
    def get_geneStrand(self):
        return self.gene_strand
                
    def add_transcript(self,sId):
        self.transcript_to_exon[sId]=[]
        
    def get_transcript(self):
        return self.transcript_to_exon.keys()
        
    def get_transcriptRelation(self,sTrId):
        return self.transcript_to_exon[sTrId]
    
    def set_transcriptRelation(self,sTrId,sExId):
        self.transcript_to_exon[sTrId].append(sExId)
    
    def update_transcriptList(self):
        for oExon in self.get_exonList():
            for sTranscript in oExon.get_transcript():
                self.set_transcriptRelation(sTranscript,oExon.get_index())
    
    def correct_exonCoord(self):
        #oFirstExon=self.get_exonList()[0]
        iAbsoluStart=self.get_geneStart()#oFirstExon.get_start()
        for oExon in self.get_exonList():
            iRelativeStart=oExon.get_start()-iAbsoluStart
            iRelativeStop=oExon.get_stop()-iAbsoluStart
            oExon.set_start(iRelativeStart)
            oExon.set_stop(iRelativeStop)
    
    def sort_exonList(self):
        tTemp=[]
        dDict={}
        tListExon=self.get_exonList()
        for oExon in tListExon:
            try:
                dDict[oExon.get_start()].append(oExon)
            except KeyError:
                dDict[oExon.get_start()]=[oExon]
        for iStart in sorted(dDict.keys()):
            for oExon in sorted(dDict[iStart], key=lambda x: x.stop):
                tTemp.append(oExon)
        self.set_exonList(tTemp)
    
    def assign_exonId(self):
        dTrId2Order={}
        dTrId2ExonNumber={}
        iExonId=0
        for oExon in self.get_exonList():
            iExonId+=1
            for iTrId in oExon.get_transcript():
                if iTrId not in dTrId2Order:
                    dTrId2Order[iTrId]=len(dTrId2Order)+1
                if iTrId not in dTrId2ExonNumber:
                    dTrId2ExonNumber[iTrId]=1
                else:
                    dTrId2ExonNumber[iTrId]+=1
            else:
                sTrName=EXON_BASENAME+str(iExonId)
                sTrName+=EXON_SEPARATOR
                sTrName+=EXON_INTERNALSEPARATOIR.join(sorted([str(dTrId2Order[X]) for X in oExon.get_transcript()]))
            oExon.set_index(sTrName)
            
    def set_emptyExonList(self):
        self.exon_list=[]
        
    def get_exonList(self):
        return self.exon_list
        
    def set_exonList(self,tList):
        self.exon_list=list(tList)
        
    def add_exon(self,oExon):
        self.exon_list.append(oExon)
        
    def remove_exon(self,oExon):
        self.exon_list.remove(oExon)
        
    def update_exonList(self,oNewExon):
        iNewStart=oNewExon.get_start()
        iNewStop=oNewExon.get_stop()
        tExonList=self.get_exonList()
        bTheSame=False
        for oOtherExon in tExonList:
            if iNewStart==oOtherExon.get_start() and iNewStop==oOtherExon.get_stop():
                oOtherExon.update(oNewExon)
                bTheSame=True
                break
        if not bTheSame:
            self.add_exon(oNewExon)
        
    def describe_exon(self):
        sContent=""
        for oExon in self.get_exonList():
            if len(sContent)!=0:
                sContent+"\n"
            sContent+=oExon.describe_self()
        return sContent

    def describe_transcript(self):
        sContent=""
        for sTrId in self.get_transcript():
            if len(sContent)!=0:
                sContent+="\n"
            sContent+="{}\t{}".format(sTrId,"\t".join(self.get_transcriptRelation(sTrId)))
        print(sContent)
        return sContent


class AlignedMatrixContent():
    
    def __init__(self,sGeneId,sReadPaf,sOutputFile):
        
        self.size=None
        
        self.currentMatrix=[]
        self.matrix_lineName=[]
        self.globalVector=[]
        self.globalBlockName=[]
        self.currentPopVector=[]
        self.globalPopVector=[]
        self.blockConfidence={}
        
        self.tr_dict_data={}
        self.parse_pafFile(sReadPaf,sGeneId)
        self.setup_matrix()
        
        self.setup_popVector()
        self.setup_globalVector()
        self.assign_blockName()
        self.assign_blockConfidence()
        #print(self.translate_cgalcodeGrammar(self.get_globalVector()))
        
        #self.print_globalVector(4)
        #self.print_globalVector(3)
        #self.print_globalVector(2)
        #self.print_globalVector(1)
        
        #iDebugLineIndex=1
        #iStartDebug=13350
        #iStopDebug=13450
        #iStartDebug=14430
        #iStopDebug=14475
        #print(self.get_matrix_lineName(iDebugLineIndex))
        #print(self.get_submatrix_byLine(iDebugLineIndex,iDebugLineIndex))
        #exit()
        
        #print(self.get_globalVector()[iStartDebug:iStopDebug+1])
        #print(self.get_globalPopVector()[iStartDebug:iStopDebug+1])
        #print(self.get_submatrix_byCol(iStartDebug,iStopDebug)[iDebugLineIndex])
        
        self.regroup_globalVector()  #This will modify the matrix!!
        
        self.setup_popVector()
        self.setup_globalVector()
        self.assign_blockName()
        self.assign_blockConfidence()
        #print(self.get_globalVector()[iStartDebug:iStopDebug+1])
        #print(self.get_globalPopVector()[iStartDebug:iStopDebug+1])
        #print(self.get_submatrix_byCol(iStartDebug,iStopDebug)[iDebugLineIndex])
        #exit()
        
        #print(self.get_line_BlockName_Structure(iDebugLineIndex))
        #exit()
        
        #print(self.get_limitedBlockName())
        #print(self.get_blockSize())
        #print(self.get_blockCoord())
        #print(self.get_blockPop())
        #for iLineIndex in range(len(self.get_matrix())):
            #print(self.get_matrix_lineName(iLineIndex))
            ##print(self.translate_cgalcodeGrammar(self.get_lineName_Vector(iLineIndex),True))
            #print(self.get_line_BlockName_Structure(iLineIndex))
            ##self.print_individualVector(self.get_lineName_Vector(iLineIndex))
            ##print(self.get_lineName_Vector(iLineIndex))
            ##exit()
            
        self.globalData2tsv(sOutputFile)
    
    def get_groupOfBlock(self):
        tBlockName=self.get_line_blockName()
        dBlock2Group={}
        tCurrentGroup=[]
        sPreviousBlockName=""
        for sBlockName in tBlockName:
            if sPreviousBlockName==sBlockName:
                continue
            if sBlockName=="i":
                if len(tCurrentGroup)==0:
                    continue
                for sSoloBlock in tCurrentGroup:
                    dBlock2Group[sSoloBlock]=tuple(tCurrentGroup)
                tCurrentGroup=[]
            else:
                tCurrentGroup.append(sBlockName)
            sPreviousBlockName=sBlockName
        return dBlock2Group
    
    def assign_blockConfidence(self):
        dBlock2Group=self.get_groupOfBlock()
        tListOfGroup=sorted(list(set(dBlock2Group.values())))
        dBlock2Pop=self.get_blockPop()
        dBlock2Confidence={}
        for tGroup in tListOfGroup:
            iMaxPop=0
            for sBlock in tGroup:
                iMaxPop=max(iMaxPop,dBlock2Pop[sBlock])
            for sBlock in tGroup:
                fConfidence=round(dBlock2Pop[sBlock]/iMaxPop*100,2)
                dBlock2Confidence[sBlock]=fConfidence
        self.blockConfidence=dBlock2Confidence
        
    def get_blockConfidence(self,sBlockName=None):
        if sBlockName is None:
            return self.blockConfidence
        return self.blockConfidence[sBlockName]
    
    def globalData2tsv(self,sOutputFile="Default_globalData2tsv.tsv"):
        tBlockName=self.get_limitedBlockName()
        dBlockName2Size=self.get_blockSize()
        dBlockName2Coord=self.get_blockCoord()
        dBlockName2Pop=self.get_blockPop()
        dBlockName2Confidence=self.get_blockConfidence()
        sHeader="ReadId\t{}\n".format("\t".join(tBlockName))
        sStartLine="Start"
        sStopLine="Stop"
        sPopLine="Population"
        sSizeLine="Size"
        sConfidenceLine="%Confidence"
        for sBlockName in tBlockName:
            sStartLine+="\t"
            sStopLine+="\t"
            sPopLine+="\t"
            sSizeLine+="\t"
            sConfidenceLine+="\t"
            if sBlockName=="i":
                continue
            sStartLine+="{}".format(dBlockName2Coord[sBlockName][0])
            sStopLine+="{}".format(dBlockName2Coord[sBlockName][-1])
            sPopLine+="{}".format(dBlockName2Pop[sBlockName])
            sSizeLine+="{}".format(dBlockName2Size[sBlockName])
            sConfidenceLine+="{}".format(dBlockName2Confidence[sBlockName])
        sStartLine+="\n"
        sStopLine+="\n"
        sPopLine+="\n"
        sSizeLine+="\n"
        sConfidenceLine+="\n"
        sCoreLine=""
        for iLineIndex in range(len(self.get_matrix())):
            tLineModel=self.get_line_BlockName_Structure(iLineIndex)
            print(self.get_matrix_lineName(iLineIndex))
            print(tLineModel)
            sCoreLine+=self.get_matrix_lineName(iLineIndex)
            tData=[]
            for sBlockName in tBlockName:
                if sBlockName=="i":
                    tData.append("")
                elif sBlockName in tLineModel:
                    tData.append("X")
                else:
                    tData.append("")
            iPreviousValue=None
            iPreviousIndex=None
            for iInLineIndex in range(len(tLineModel)):
                if isinstance(tLineModel[iInLineIndex],int):
                    if iInLineIndex==0:
                        for iThisIndex in range(len(tData)):
                            if tData[iThisIndex]=="X":
                                tData[iThisIndex-1]=str(tLineModel[iInLineIndex])
                                break
                    elif iInLineIndex==len(tLineModel)-1 or (iInLineIndex==len(tLineModel)-2 and tLineModel[-1]=="i"):
                        for iReverseIndex in range(len(tData)-1,-1,-1):
                            if tData[iReverseIndex]=="X":
                                tData[iReverseIndex+1]=str(tLineModel[iInLineIndex])
                                break
                    else:
                        if iPreviousValue is None:
                            iPreviousValue=tLineModel[iInLineIndex]
                            iPreviousIndex=iInLineIndex
                        elif iPreviousValue==tLineModel[iInLineIndex]:
                            iStartingMissingIndex=tBlockName.index(tLineModel[iPreviousIndex-1])+1
                            iEndingMissingIndex=tBlockName.index(tLineModel[iInLineIndex+1])-1
                            iValueMissingIndex=int((iStartingMissingIndex+iEndingMissingIndex)/2)
                            for iAnotherIndex in range(iStartingMissingIndex,iEndingMissingIndex+1):
                                tData[iAnotherIndex]="****"
                                if iAnotherIndex==iValueMissingIndex:
                                    tData[iAnotherIndex]=str(iPreviousValue)
                            iPreviousValue=None
                            iPreviousIndex=None
                        else:
                            print(iPreviousValue,tLineModel[iInLineIndex])
                            exit("Error 1805 : non-equivalent missing part")
                            
                        
                        
                    
            sCoreLine+="\t"+"\t".join(list(tData))+"\n"
        FILE=open(sOutputFile,"w")
        FILE.write(sStartLine+sStopLine+sSizeLine+sHeader+sCoreLine+sPopLine+sConfidenceLine)
        FILE.close()
        
    
    def get_1Dmatrix(self):
        return [bool(X) for X in self.get_globalPopVector()]
        
    def get_blockSize(self,sBlockName=None):
        tBlockName=self.get_globalBlockName()
        if sBlockName is not None:
            return len([X for X in tBlockName if X==sBlockName])
        dBlock2Size={}
        for sName in tBlockName:
            if sName=="i":
                continue
            try:
                dBlock2Size[sName]+=1
            except KeyError:
                dBlock2Size[sName]=1
        return dBlock2Size
    
    def get_blockPop(self,sBlockName=None):
        tGlobalPop=self.get_globalPopVector()
        tBlockName=self.get_globalBlockName()
        dResult={}
        for iIndex in range(len(tBlockName)):
            dResult[tBlockName[iIndex]]=tGlobalPop[iIndex]
        if sBlockName is None:
            return dResult
        return dResult[sBlockName]
    
    def get_blockCoord(self,sBlockName=None):
        tBlockName=self.get_globalBlockName()
        dBlock2Coord={}
        iStartValue=-1
        iStopValue=-1
        sPreviousName=None
        for iIndex in range(len(tBlockName)):
            sName=tBlockName[iIndex]
            if sName=="i":
                continue
            if sPreviousName is None:
                iStartValue=iIndex
                iStopValue=iIndex
                sPreviousName=sName
            else:
                if sName==sPreviousName:
                    iStopValue+=1
                else:
                    dBlock2Coord[sPreviousName]=(iStartValue,iStopValue)
                    iStartValue=iIndex
                    iStopValue=iIndex
                    sPreviousName=sName
        if sPreviousName not in dBlock2Coord:
            dBlock2Coord[sPreviousName]=(iStartValue,iStopValue)
        
        if sBlockName is None:
            return dBlock2Coord
        return dBlock2Coord[sBlockName]
    
    def get_line_blockName(self,iLineIndex=None):
        tIndividualBlockNameVector=[]
        if iLineIndex is not None:
            tTargetedMatrixLine=self.get_matrix()[iLineIndex]
        else:
            tTargetedMatrixLine=self.get_1Dmatrix()
        for iColIndex in range(len(tTargetedMatrixLine)):
            if tTargetedMatrixLine[iColIndex]==1:
                tIndividualBlockNameVector.append(self.get_globalBlockName(iColIndex))
            elif tTargetedMatrixLine[iColIndex]==0:
                tIndividualBlockNameVector.append("i")
            elif tTargetedMatrixLine[iColIndex]<0:
                tIndividualBlockNameVector.append(tTargetedMatrixLine[iColIndex])
        return tIndividualBlockNameVector
    
    def get_line_BlockName_Structure(self,iLineIndex=None):
        #for iColIndex in range(len(self.get_globalBlockName())):
            #print(iColIndex,self.get_globalVector(iColIndex),self.get_globalBlockName(iColIndex),self.get_matrix()[iLineIndex][iColIndex])
        tIndividualBlockNameVector=self.get_line_blockName(iLineIndex)
        #print(tIndividualBlockNameVector)
        tPrintBlockName=[]
        for iColIndex in range(len(tIndividualBlockNameVector)):
            #print(tIndividualBlockNameVector[iColIndex],)
            if iColIndex==0:
                if tIndividualBlockNameVector[iColIndex]=="i":
                    continue
                else:
                    tPrintBlockName.append(tIndividualBlockNameVector[iColIndex])
            elif tIndividualBlockNameVector[iColIndex]==tIndividualBlockNameVector[iColIndex-1]:
                continue
            else:
                tPrintBlockName.append(tIndividualBlockNameVector[iColIndex])
        #print(tPrintBlockName)
        return tPrintBlockName
            
    def get_limitedBlockName(self):
        tBlockLine=self.get_line_blockName()
        tLimitedBlock=[]
        for iIndex in range(len(tBlockLine)):
            if iIndex==0:
                tLimitedBlock.append(tBlockLine[iIndex])
                continue
            if tBlockLine[iIndex]!=tBlockLine[iIndex-1]:
                tLimitedBlock.append(tBlockLine[iIndex])
        return tLimitedBlock
        
    
    def get_lineName_Vector(self,iLineIndex):
        tIndividualVector=[]
        tGlobalVector=self.get_globalVector()
        tTargetedMatrixLine=self.get_matrix()[iLineIndex]
        for iColIndex in range(len(tGlobalVector)):
            if tTargetedMatrixLine[iColIndex]==1:
                tIndividualVector.append(tGlobalVector[iColIndex])
                #print(tTargetedMatrixLine[iColIndex],tGlobalVector[iColIndex])
            else:
                tIndividualVector.append("i")
                #print(tTargetedMatrixLine[iColIndex],"Nothing")
        #print(tIndividualVector)
        return tIndividualVector
            
    def assign_blockName(self):
        tGlobalVector=self.get_globalVector()
        tBlockNameVector=list(tGlobalVector)
        sBlockName=None
        for iIndex in range(len(tGlobalVector)):
            if tGlobalVector[iIndex] not in ["i","="]:
                sBlockName=self.update_blockName(sBlockName)
                tBlockNameVector[iIndex]=sBlockName
            elif tGlobalVector[iIndex]=="=":
                tBlockNameVector[iIndex]=sBlockName
            #print(iIndex,tGlobalVector[iIndex],tBlockNameVector[iIndex])
        self.set_globalBlockName(list(tBlockNameVector))
        
    def set_globalBlockName(self,oValue,iIndex=None):        
        if iIndex is None:
            self.globalBlockName=list(oValue)
        else:
            self.globalBlockName[iIndex]=oValue
            
    def get_globalBlockName(self,iIndex=None):
        if iIndex is None:
            return self.globalBlockName
        return self.globalBlockName[iIndex]
        
    def update_blockName(self,sName=None):
        if sName is None:
            return "A"
        if sName=="Z":
            return "A1"
        if len(sName)!=1:
            return sName[0]+str((int(sName[1:])+1))
        return chr(ord(sName)+1)
    
    def regroup_globalVector(self):
        #REGROUP_GLOBALVECTOR_VALUE
        iLastIndex=None
        iCurrentGroup=0
        dGroup2Index={}
        for iIndex in range(len(self.get_globalVector())):
            iAddIndex=None
            if self.get_globalVector(iIndex) not in ["i","="]:
                iAddIndex=iIndex
            if self.get_globalVector(iIndex)=="i" and self.get_globalVector(iIndex-1)=="=":
                 iAddIndex=iIndex-1
            #if self.get_globalVector(iIndex) not in ["i","="] or (self.get_globalVector(iIndex)=="i" and self.get_globalVector(iIndex-1)=="="):
            if iAddIndex is not None:
                if iLastIndex is None:
                    iLastIndex=iAddIndex
                    continue
                if iLastIndex+REGROUP_GLOBALVECTOR_VALUE>=iAddIndex:
                    #print(iAddIndex,"grouped")
                    if iCurrentGroup not in dGroup2Index:
                        dGroup2Index[iCurrentGroup]=[iLastIndex,iAddIndex]
                    else:
                        dGroup2Index[iCurrentGroup].append(iAddIndex)
                else:
                    #print(iAddIndex,"not grouped")
                    iCurrentGroup+=1
                iLastIndex=iAddIndex
        #print(dGroup2Index)
        
        ##ASSIGN_GROUP_VALUE
        for iGroupId in sorted(dGroup2Index):
            print("iGroupId",dGroup2Index[iGroupId])
            if len(set([self.get_globalVector(X) for X in dGroup2Index[iGroupId] if self.get_globalVector(X) in [">","<"]]))!=1:
                #exit("ERROR 1726 : different type of variation\n{}".format([self.get_globalVector(X) for X in dGroup2Index[iGroupId] if self.get_globalVector(X) in [">","<"]]))
                print("NEAR-ERROR 1726 : different type of variation\n{}\t{}".format([self.get_globalVector(X) for X in dGroup2Index[iGroupId] if self.get_globalVector(X) in [">","<"]],dGroup2Index[iGroupId]))
                print("Try without group this elements")
                continue
            else:
                sTag=[self.get_globalVector(X) for X in dGroup2Index[iGroupId] if self.get_globalVector(X) in [">","<"]][0]
            if len([self.get_globalVector(X) for X in dGroup2Index[iGroupId] if self.get_globalVector(X) not in [">","<"]])==1:
                sBlockId=[self.get_globalVector(X) for X in dGroup2Index[iGroupId] if self.get_globalVector(X) not in [">","<"]][0]
            else:
                sBlockId=sTag
            print("sBlockId",sBlockId,"sTag",sTag)
            iMajorityIndex=None
            iMajorityPop=None
            tEqualityIndex=[]
            iSumPop=0
            for iIndex in dGroup2Index[iGroupId]:
                if iIndex==0:
                    print("NEAR-ERROR 1742 : iIndex=0, no iIndex-1 possible")
                    print("Maybe this will work...")
                if sTag=="<":
                    iPop=self.get_deltaPop_forNewStart(iIndex)
                elif sTag==">":
                    iPop=self.get_deltaPop_forNewStop(iIndex)
                else:
                    #Starting block, same as <
                    iPop=self.get_currentPopVector(iIndex)-self.get_currentPopVector(iIndex-1)
                iSumPop+=iPop
                #print(iIndex,self.get_globalVector(iIndex),self.get_globalPopVector(iIndex),iPop)
                #print(range(iIndex-5,iIndex+5))
                #print(self.get_globalPopVector()[iIndex-5:iIndex+5])
                #print(self.get_globalVector()[iIndex-5:iIndex+5])
                
                if iMajorityIndex is None:
                    iMajorityIndex=iIndex
                    iMajorityPop=iPop
                    tEqualityIndex=[iIndex]
                elif iPop>iMajorityPop:
                    iMajorityIndex=iIndex
                    iMajorityPop=iPop
                    tEqualityIndex=[iIndex]
                elif iPop==iMajorityPop:
                    tEqualityIndex.append(iIndex)
            #print(iGroupId,iMajorityIndex)
            
            
            if len(tEqualityIndex)!=1:    
                #print("MajIndex",iMajorityIndex,"Pop",iPop,"Equality",tEqualityIndex)
                #print("-------------")
                if sTag=="<":
                    #cumulative start, take the last
                    iMajorityIndex=tEqualityIndex[-1]
                else:
                    #cumulative end, take the first
                    iMajorityIndex=tEqualityIndex[0]
                print("MajIndex",iMajorityIndex,"Pop",iPop,"Equality",tEqualityIndex)
            
            ##DEBUG
            #print("Index",dGroup2Index[iGroupId])
            #if sTag=="<":
                #print("realPop",[self.get_deltaPop_forNewStart(X) for X in dGroup2Index[iGroupId]])
            #else:
                #print("realPop",[self.get_deltaPop_forNewStop(X) for X in dGroup2Index[iGroupId]])
            #print("MajIndex",iMajorityIndex,"majPop",iMajorityPop,"allPop",iSumPop)
            #print("globalPop",[self.get_globalPopVector(X) for X in dGroup2Index[iGroupId]])
            #print("nearGlobalPop",[self.get_globalPopVector(X) for X in [dGroup2Index[iGroupId][0]-1,dGroup2Index[iGroupId][-1]+1]])
            #print("globalVector",[self.get_globalVector(X) for X in dGroup2Index[iGroupId]])
            #print("nearGlobalVector",[self.get_globalVector(X) for X in [dGroup2Index[iGroupId][0]-1,dGroup2Index[iGroupId][-1]+1]])
            #print("-------------")
            
            ##MERGE_GROUP
            #print(dGroup2Index[iGroupId],"become",sTag,iMajorityIndex,iSumPop)
            #print("before",self.get_submatrix_byCol(dGroup2Index[iGroupId][0],dGroup2Index[iGroupId][-1]))
            tMatrix=self.get_matrix()
            iUpstreamColIndex=dGroup2Index[iGroupId][0]-1
            iDownstreamColIndex=dGroup2Index[iGroupId][-1]+1
            iStartColIndex=dGroup2Index[iGroupId][0]
            iStopColIndex=dGroup2Index[iGroupId][-1]
            if sTag=="<":
                iStartColIndex=dGroup2Index[iGroupId][0]
                iStopColIndex=dGroup2Index[iGroupId][-1]
                for iLineIndex in range(len(tMatrix)):
                    #print(self.get_matrix_lineName()[iLineIndex])
                    if tMatrix[iLineIndex][iUpstreamColIndex]==1:
                        #This read is already aligned before this section and is not concerned by changement
                        continue
                    if len([X for X in tMatrix[iLineIndex][iStartColIndex:iStopColIndex+1] if X==1])==0:
                        #This read is not aligned during this section and is not concerned by changement
                        continue
                    #print("A>>>",tMatrix[iLineIndex][iStartColIndex-2:iStopColIndex+2])
                    bHaveNegativeValue=False
                    iNegativeValue=None
                    iNegativeValueIndex=None
                    #print("B>>>",tMatrix[iLineIndex][iStartColIndex-2:iStopColIndex+2])
                    for iColIndex in range(iStartColIndex-1,iStopColIndex+1): # from iStartColIndex-1 to take into account -40[< 1 1 1 1]
                        if tMatrix[iLineIndex][iColIndex]<0:
                            bHaveNegativeValue=True
                            iNegativeValue=tMatrix[iLineIndex][iColIndex]
                            iNegativeValueIndex=iColIndex
                            tMatrix[iLineIndex][iColIndex]=0 # turn it to 0 now
                            break
                            #print(">",iNegativeValue,iNegativeValueIndex)
                    #print("C>>>",tMatrix[iLineIndex][iStartColIndex-2:iStopColIndex+2])
                    for iColIndex in range(iStartColIndex,iStopColIndex+1):
                        if iColIndex<iMajorityIndex:
                            tMatrix[iLineIndex][iColIndex]=0                                
                        if iColIndex>=iMajorityIndex:
                            tMatrix[iLineIndex][iColIndex]=1
                        if bHaveNegativeValue:
                            #print(iColIndex,iMajorityIndex-1)
                            if iColIndex==iMajorityIndex-1:
                                #print("Z",tMatrix[iLineIndex][iStartColIndex:iStopColIndex+2])
                                tMatrix[iLineIndex][iColIndex]=iNegativeValue
                                #print("Z",tMatrix[iLineIndex][iStartColIndex:iStopColIndex+2])
                    #print("D>>>",tMatrix[iLineIndex][iStartColIndex-2:iStopColIndex+2])
                            
            else:
                iStartColIndex=dGroup2Index[iGroupId][0]-1
                iStopColIndex=dGroup2Index[iGroupId][-1]
                for iLineIndex in range(len(tMatrix)):
                    #print(self.get_matrix_lineName()[iLineIndex])
                    if tMatrix[iLineIndex][iDownstreamColIndex]==1:
                        #This read is ever aligned after this section and is not concerned by changement
                        #print("Apres",tMatrix[iLineIndex][iStartColIndex:iStopColIndex+2])
                        continue
                    if len([X for X in tMatrix[iLineIndex][iStartColIndex:iStopColIndex+1] if X==1])==0:
                        #This read is not aligned during this section and is not concerned by changement
                        #print("Apres",tMatrix[iLineIndex][iStartColIndex:iStopColIndex+2])
                        continue
                    #print("A'>>>",tMatrix[iLineIndex][iStartColIndex-2:iStopColIndex+2])
                    bHaveNegativeValue=False
                    iNegativeValue=None
                    iNegativeValueIndex=None
                    #print("B'>>>",tMatrix[iLineIndex][iStartColIndex-2:iStopColIndex+2])
                    for iColIndex in range(iStartColIndex,iStopColIndex+1):
                        if tMatrix[iLineIndex][iColIndex]<0:
                            bHaveNegativeValue=True
                            iNegativeValue=tMatrix[iLineIndex][iColIndex]
                            iNegativeValueIndex=iColIndex
                    #print("C'>>>",tMatrix[iLineIndex][iStartColIndex-2:iStopColIndex+2])
                    for iColIndex in range(iStartColIndex,iStopColIndex+1):
                        if iColIndex<=iMajorityIndex:
                            tMatrix[iLineIndex][iColIndex]=1
                        if iColIndex>iMajorityIndex:
                            tMatrix[iLineIndex][iColIndex]=0
                        if bHaveNegativeValue:
                            if iColIndex==iMajorityIndex+1:
                                tMatrix[iLineIndex][iColIndex]=iNegativeValue
                    #print("D'>>>",tMatrix[iLineIndex][iStartColIndex-2:iStopColIndex+2])
                    
            if sTag=="<":
                sPreviousTag=self.get_globalVector(dGroup2Index[iGroupId][0]-1)
                iPreviousPop=self.get_globalPopVector(dGroup2Index[iGroupId][0]-1)
                for iIndex in range(dGroup2Index[iGroupId][0],dGroup2Index[iGroupId][-1]+1):
                    if iIndex<iMajorityIndex:
                        self.set_globalVector(sPreviousTag,iIndex)
                        self.set_globalPopVector(iPreviousPop,iIndex)
                    if iIndex==iMajorityIndex:
                        self.set_globalVector(sBlockId,iIndex)
                        self.set_globalPopVector(iSumPop,iIndex)
                    if iIndex>iMajorityIndex:
                        self.set_globalVector("=",iIndex)
                        self.set_globalPopVector(iSumPop,iIndex)
                    #print(self.get_globalVector(iIndex),self.get_globalPopVector(iIndex))
            else:
                sPreviousTag=self.get_globalVector(dGroup2Index[iGroupId][0]-1)
                iPreviousPop=self.get_globalPopVector(dGroup2Index[iGroupId][0]-1)
                sNextTag=self.get_globalVector(dGroup2Index[iGroupId][0]+1)
                #print("PreviousTag/Pop",sPreviousTag,iPreviousPop)
                for iIndex in range(dGroup2Index[iGroupId][0],dGroup2Index[iGroupId][-1]+1):
                    if iIndex<iMajorityIndex:
                        self.set_globalVector(sPreviousTag,iIndex)
                        self.set_globalPopVector(iPreviousPop,iIndex)
                    if iIndex==iMajorityIndex:
                        self.set_globalVector(sBlockId,iIndex)
                        self.set_globalPopVector(iSumPop,iIndex)
                    if iIndex>iMajorityIndex:
                        if sNextTag!="i":
                            self.set_globalVector("=",iIndex)
                            self.set_globalPopVector(iSumPop,iIndex)
                        else:
                            self.set_globalVector("i",iIndex)
                            self.set_globalPopVector(0,iIndex)
                    #print(self.get_globalVector(iIndex),self.get_globalPopVector(iIndex))
            #print("after",self.get_submatrix_byCol(dGroup2Index[iGroupId][0],dGroup2Index[iGroupId][-1]))
            #print("Change")
            #print("MajIndex",iMajorityIndex,"majPop",iMajorityPop,"allPop",iSumPop)
            #print("globalPop",[self.get_globalPopVector(X) for X in dGroup2Index[iGroupId]])
            #print("nearGlobalPop",[self.get_globalPopVector(X) for X in [dGroup2Index[iGroupId][0]-1,dGroup2Index[iGroupId][-1]+1]])
            #print("globalVector",[self.get_globalVector(X) for X in dGroup2Index[iGroupId]])
            #print("nearGlobalVector",[self.get_globalVector(X) for X in [dGroup2Index[iGroupId][0]-1,dGroup2Index[iGroupId][-1]+1]])
            #print("-------------")
    
    def print_individualVector(self,tVector):
        tList=[]
        sPreviousChar=None
        for iIndex in range(len(tVector)):
            if iIndex==0:
                tList.append(tVector[iIndex])
            elif iIndex==len(tVector)-1:
                tList.append(tVector[iIndex])
            elif tVector[iIndex]=="=":
                continue
            elif sPreviousChar!=tVector[iIndex]:
                tList.append(tVector[iIndex])
            elif sPreviousChar==tVector[iIndex] and sPreviousChar in [">","<"] and tVector[iIndex] in [">","<"]:
                tList.append(tVector[iIndex])
            sPreviousChar=tVector[iIndex]
        print(tList)
    
    def print_globalVector(self,iDebugMode=0):
        tGlobalVector=self.get_globalVector()
        if iDebugMode==0:
            print(tGlobalVector)
        elif iDebugMode==1:
            for iIndex in range(len(tGlobalVector)):
                #print("L{}".format(iIndex+1),tGlobalVector[iIndex])
                print("L{}".format(iIndex+1),tGlobalVector[iIndex],self.get_currentPopVector(iIndex),
                        "p{}".format(self.get_deltaPop_wPrevious(iIndex)),"n{}".format(self.get_deltaPop_wNext(iIndex)))
        elif iDebugMode==2:
            sPreviousChar=None
            for iIndex in range(len(tGlobalVector)):
                if iIndex==0:
                    print("L{}".format(iIndex+1),tGlobalVector[iIndex],self.get_currentPopVector(iIndex),
                            "p{}".format(self.get_deltaPop_wPrevious(iIndex)),"n{}".format(self.get_deltaPop_wNext(iIndex)))
                elif iIndex==len(tGlobalVector)-1:
                    print("L{}".format(iIndex+1),tGlobalVector[iIndex],self.get_currentPopVector(iIndex),
                            "p{}".format(self.get_deltaPop_wPrevious(iIndex)),"n{}".format(self.get_deltaPop_wNext(iIndex)))
                elif sPreviousChar!=tGlobalVector[iIndex]:
                    print("L{}".format(iIndex+1),tGlobalVector[iIndex],self.get_currentPopVector(iIndex),
                            "p{}".format(self.get_deltaPop_wPrevious(iIndex)),"n{}".format(self.get_deltaPop_wNext(iIndex)))
                elif sPreviousChar==tGlobalVector[iIndex] and sPreviousChar in [">","<"] and tGlobalVector[iIndex] in [">","<"]:
                    print("L{}".format(iIndex+1),tGlobalVector[iIndex],self.get_currentPopVector(iIndex),
                            "p{}".format(self.get_deltaPop_wPrevious(iIndex)),"n{}".format(self.get_deltaPop_wNext(iIndex)))
                sPreviousChar=tGlobalVector[iIndex]
        elif iDebugMode==3:
            tList=[]
            tListPop=[]
            sPreviousChar=None
            for iIndex in range(len(tGlobalVector)):
                if iIndex==0:
                    tList.append(tGlobalVector[iIndex])
                    tListPop.append(self.get_currentPopVector(iIndex))
                elif iIndex==len(tGlobalVector)-1:
                    tList.append(tGlobalVector[iIndex])
                    tListPop.append(self.get_currentPopVector(iIndex))
                elif sPreviousChar!=tGlobalVector[iIndex]:
                    tList.append(tGlobalVector[iIndex])
                    tListPop.append(self.get_currentPopVector(iIndex))
                elif sPreviousChar==tGlobalVector[iIndex] and sPreviousChar in [">","<"] and tGlobalVector[iIndex] in [">","<"]:
                    tList.append(tGlobalVector[iIndex])
                    tListPop.append(self.get_currentPopVector(iIndex))
                sPreviousChar=tGlobalVector[iIndex]
            print(tList)
            print(tListPop)
        elif iDebugMode==4:
            tList=[]
            tListIndex=[]
            sPreviousChar=None
            for iIndex in range(len(tGlobalVector)):
                if iIndex==0:
                    tList.append(tGlobalVector[iIndex])
                    tListIndex.append(iIndex)
                elif iIndex==len(tGlobalVector)-1:
                    tList.append(tGlobalVector[iIndex])
                    tListIndex.append(iIndex)
                elif sPreviousChar!=tGlobalVector[iIndex]:
                    tList.append(tGlobalVector[iIndex])
                    tListIndex.append(iIndex)
                elif sPreviousChar==tGlobalVector[iIndex] and sPreviousChar in [">","<"] and tGlobalVector[iIndex] in [">","<"]:
                    tList.append(tGlobalVector[iIndex])
                    tListIndex.append(iIndex)
                sPreviousChar=tGlobalVector[iIndex]
            print(tList)
            print(tListIndex)
        elif iDebugMode==5:
            tList=[]
            tListIndex=[]
            sPreviousChar=None
            for iIndex in range(len(tGlobalVector)):
                if iIndex==0:
                    tList.append(tGlobalVector[iIndex])
                    tListIndex.append(iIndex)
                elif iIndex==len(tGlobalVector)-1:
                    tList.append(tGlobalVector[iIndex])
                    tListIndex.append(iIndex)
                elif tGlobalVector[iIndex]=="=":
                    continue
                elif sPreviousChar!=tGlobalVector[iIndex]:
                    tList.append(tGlobalVector[iIndex])
                    tListIndex.append(iIndex)
                elif sPreviousChar==tGlobalVector[iIndex] and sPreviousChar in [">","<"] and tGlobalVector[iIndex] in [">","<"]:
                    tList.append(tGlobalVector[iIndex])
                    tListIndex.append(iIndex)
                sPreviousChar=tGlobalVector[iIndex]
            print(tList)
            print(tListIndex)

    def setup_popVector(self):
        tMatrix=self.get_matrix()
        tGlobalPopVector=[]
        for iColIndex in range(len(tMatrix[0])):
            iPop=0
            for iLineIndex in range(len(tMatrix)):
                if tMatrix[iLineIndex][iColIndex]>=0:
                    iPop+=tMatrix[iLineIndex][iColIndex]
                else:
                    iPop+=0 #-X value for unaligned nt are assumed as 0 pop
            tGlobalPopVector.append(iPop)
        self.globalPopVector=list(tGlobalPopVector)
    
    def get_deltaPop_forNewStart(self,iIndex):
        if iIndex!=0:
            iPreviousPop=self.get_globalPopVector(iIndex-1)
            iCurrentPop=self.get_globalPopVector(iIndex)
            return iCurrentPop-iPreviousPop
        else:
            return 0
    
    def get_deltaPop_forNewStop(self,iIndex):
        if iIndex!=len(self.get_globalPopVector())-1:
            iPreviousPop=self.get_globalPopVector(iIndex-1)
            iCurrentPop=self.get_globalPopVector(iIndex)
            if iPreviousPop==iCurrentPop: 
                return iCurrentPop
            return iPreviousPop-iCurrentPop
        else:
            return 0
        #if iIndex!=len(self.get_globalPopVector())-1:
            #iCurrentPop=self.get_globalPopVector(iIndex)
            #for iNextIndex in range(iIndex,len(self.get_globalPopVector())):
                #if self.get_globalPopVector(iNextIndex)!=iCurrentPop:
                    #break
            #iNextPop=self.get_globalPopVector(iNextIndex)
            #return iCurrentPop-iNextPop
        #else:
            #return 0
    
    def get_globalPopVector(self,iIndex=None):
        if iIndex is None:
            return self.globalPopVector
        else:
            return self.globalPopVector[iIndex]
    
    def set_globalPopVector(self,oValue,iIndex=None):
        if iIndex is None:
            self.globalPopVector=list(oValue)
        else:
            self.globalPopVector[iIndex]=oValue
    
    def setup_globalVector(self):
        iSize=self.get_size()
        tVector=['']*iSize
        tMatrix=self.get_matrix()
        iBlockIndex=0
        dbPreviousVector=None
        for iColIndex in range(len(tVector)):
            tCurrentVector=[]
            for iLineIndex in range(len(tMatrix)):
                if tMatrix[iLineIndex][iColIndex]>=0:
                    tCurrentVector.append(tMatrix[iLineIndex][iColIndex])
                else:
                    tCurrentVector.append(0) #-X value for unaligned nt are assumed as 0 pop
            
            dbCurrentVector=tuple(tCurrentVector)
            
            if dbPreviousVector is None:
                if sum(dbCurrentVector)==0:
                    tVector[iColIndex]="i"
                else:
                    iBlockIndex+=1
                    tVector[iColIndex]=str(iBlockIndex)
                dbPreviousVector=dbCurrentVector
                continue
                
            bIsSimilar=dbPreviousVector==dbCurrentVector
            
            bPreviousIsIncludeIntoCurrent=None
            tPreviousIsIncludeIntoCurrent=[dbPreviousVector[X] for X in range(len(dbPreviousVector)) if dbPreviousVector[X]==1 and dbPreviousVector[X]!=dbCurrentVector[X]]
            if len(tPreviousIsIncludeIntoCurrent)==0:
                bPreviousIsIncludeIntoCurrent=False
            else:
                bPreviousIsIncludeIntoCurrent=True
                
            bCurrentIsIncludeIntoPrevious=None
            tCurrentIsIncludeIntoPrevious=[dbCurrentVector[X] for X in range(len(dbCurrentVector)) if dbCurrentVector[X]==1 and dbCurrentVector[X]!=dbPreviousVector[X]]
            if len(tCurrentIsIncludeIntoPrevious)==0:
                bCurrentIsIncludeIntoPrevious=False
            else:
                bCurrentIsIncludeIntoPrevious=True

            bIsIntron=False
            if sum(dbCurrentVector)==0:
                bIsIntron=True
                
            bPreviousIsIntron=False
            if sum(dbPreviousVector)==0:
                bPreviousIsIntron=True

            sMark=""
            if bIsIntron:
                sMark="i"
            elif bIsSimilar:
                sMark="="
            elif bCurrentIsIncludeIntoPrevious and not bPreviousIsIntron:
                sMark="<"
            elif bPreviousIsIncludeIntoCurrent and not bPreviousIsIntron:
                sMark=">"
            else:
                iBlockIndex+=1
                sMark=str(iBlockIndex)
            tVector[iColIndex]=sMark
            
            dbPreviousVector=dbCurrentVector            
            
        self.globalVector=list(tVector)
                
    def get_globalVector(self,iIndex=None):
        if iIndex is None:
            return self.globalVector
        else:
            return self.globalVector[iIndex]
        
    def set_globalVector(self,oValue,iIndex=None):
        if iIndex is None:
            self.globalVector=list(oValue)
        else:
            self.globalVector[iIndex]=oValue
    
    def setup_matrix(self):
        iSize=self.get_size()
        dTr2Data=self.get_tr_dict_data()
        iLine=len(dTr2Data)
        tMatrix=[]
        tBaseVector=[0]*iSize
        for sTrId in sorted(dTr2Data):
            #print(sTrId)
            tCurrentVector=list(tBaseVector)
            self.set_matrix_lineName(sTrId,True)
            iPreviousReadStop=None
            iPreviousGeneStop=None
            iPreviousReadStart=None
            iPreviousGeneStart=None
            iLastIndex=None
            for iDbIndex in range(len(dTr2Data[sTrId])):
                dbKey=sorted(dTr2Data[sTrId])[iDbIndex]
                #print(dbKey)
                iGeneStart=dTr2Data[sTrId][dbKey]["GeneStart"]
                iGeneStop=dTr2Data[sTrId][dbKey]["GeneStop"]
                iReadStart=dTr2Data[sTrId][dbKey]["ReadStart"]
                iReadStop=dTr2Data[sTrId][dbKey]["ReadStop"]
                sStrand=dTr2Data[sTrId][dbKey]["Strand"]
                if sStrand=="+":
                    sStart="start"
                    sStop="end"
                else:
                    sStart="end"
                    sStop="start"
                #print(iGeneStart,iGeneStop,iReadStart,iReadStop,sStrand)
                if iGeneStart>iGeneStop:
                    exit("ERROR L1705 : iStart>iStop")
                for iIndex in range(iGeneStart,iGeneStop+1):
                    if iIndex<iSize:
                        tCurrentVector[iIndex]=1
                    else:
                        exit("ERROR L2257 : iIndex>=iSize")
                if iDbIndex==0:
                    if iReadStart!=0:
                        iDelta=-(iReadStart-1)
                        print("Warning : {}\n\tSoft gap at the {} of the gene ({}-{}: {}nt gap)".format(sTrId,sStart,1,iReadStart,iDelta))
                        #print("'-> corresponding to {} {} on the gene".format(iGeneStart,iGeneStop))
                        if sStrand=="+":
                            try:
                                tCurrentVector[iGeneStart-1]=iDelta
                            except IndexError:
                                tCurrentVector[iGeneStart]=iDelta+1
                        else:
                            try:
                                tCurrentVector[iGeneStop+1]=iDelta
                            except IndexError:
                                tCurrentVector[iGeneStop]=iDelta-1
                elif iDbIndex==len(dTr2Data[sTrId])-1:
                    iReadSize=dTr2Data[sTrId][dbKey]["ReadSize"]
                    if iReadStop!=iReadSize:
                        iDelta=-(iReadSize-iReadStop)
                        print("Warning : {}\n\tSoft gap at the {} of the gene ({}-{}: {}nt gap)".format(sTrId,sStop,iReadStop,iReadSize,iDelta))
                        #print("'-> corresponding to {} {} on the gene".format(iGeneStart,iGeneStop))
                        if sStrand=="+":
                            try:
                                tCurrentVector[iGeneStop+1]=iDelta
                            except IndexError:
                                tCurrentVector[iGeneStop]=iDelta+1
                        else:
                            try:
                                tCurrentVector[iGeneStart-1]=iDelta
                            except IndexError:
                                tCurrentVector[iGeneStart]=iDelta-1
                else:
                    if iPreviousReadStop+1<iReadStart:
                        iDelta=-(iReadStart-iPreviousReadStop)
                        print("Warning : {}\n\tdiscontinue alignment ({}-{}: {}nt gap)".format(sTrId,iPreviousReadStop,iReadStart,iDelta))
                        #print("'-> corresponding to {} {} on the gene".format(iGeneStart,iGeneStop))
                        if sStrand=="+":
                            tCurrentVector[iPreviousGeneStop+1]=iDelta
                            tCurrentVector[iGeneStart-1]=iDelta
                        else:
                            tCurrentVector[iPreviousGeneStart+1]=iDelta
                            tCurrentVector[iGeneStop-1]=iDelta
                iPreviousReadStop=iReadStop
                iPreviousGeneStop=iGeneStop
                iPreviousReadStart=iReadStart
                iPreviousGeneStart=iGeneStart
                iLastIndex=iIndex
            #print(tCurrentVector)
            tMatrix.append(tCurrentVector)
        self.currentMatrix=list(tMatrix)
    
    def get_matrix_lineName(self,iIndex=None):
        if iIndex is None:
            return self.matrix_lineName
        else:
            return self.matrix_lineName[iIndex]
    
    def set_matrix_lineName(self,oValue,bAdding=False):
        if bAdding:
            self.matrix_lineName.append(oValue)
        else:
            self.matrix_lineName=oValue
        
    def get_matrix(self):
        return self.currentMatrix
    
    def get_submatrix_byCol(self,iStartCol,iStopCol):
        tSubMatrix=[]
        tMatrix=self.get_matrix()
        for tLine in tMatrix:
            tSubLine=[]
            for iColIndex in range(iStartCol,iStopCol+1):
                tSubLine.append(tLine[iColIndex])
            tSubMatrix.append(list(tSubLine))
        return tSubMatrix
        
    def get_submatrix_byLine(self,iStartLine,iStopLine):
        tSubMatrix=[]
        tMatrix=self.get_matrix()
        for iLineIndex in range(iStartLine,iStopLine+1):
            tSubMatrix.append(list(tMatrix[iLineIndex]))
        return tSubMatrix
    
    def parse_pafFile(self,sPathFile,sGeneId):
        dRead2Coord={}
        for sLine in open(sPathFile):
            sLine=sLine.strip()
            tLine=sLine.split()
            #tRef=["ReadId","ReadSize","ReadStart","ReadEnd","Strand","GeneId","GeneSize","GeneStart","GeneStop","NbrNotIndel","AlignSize","Quality"]
            dRef2Content={}
            for iIndex in range(len(PAF_COLUMN)):
                try:
                    dRef2Content[PAF_COLUMN[iIndex]]=int(tLine[iIndex])
                except ValueError:
                    dRef2Content[PAF_COLUMN[iIndex]]=tLine[iIndex]
            sCurrentTrId=dRef2Content["ReadId"]
            if sCurrentTrId not in dRead2Coord:
                dRead2Coord[sCurrentTrId]={}
            dRead2Coord[sCurrentTrId][(dRef2Content["ReadStart"],dRef2Content["ReadStop"])]=dRef2Content
        self.set_size(int(dRef2Content["GeneSize"]))
        self.set_tr_dict_data(dRead2Coord)

    def set_tr_dict_data(self,dDict):
        self.tr_dict_data=dDict
    
    def get_tr_dict_data(self):
        return self.tr_dict_data

    def get_size(self):
        return self.size
        
    def set_size(self,iInt):
        self.size=iInt

    def set_geneCoord(self,iStart,iEnd,iStrand):
        self.gene_start=iStart
        self.gene_end=iEnd
        self.gene_strand=iStrand
        

class OldReadContent(TranscriptContent):
    
    def __init__(self,sGeneId,sRefFile): #,bRegroupByValue=REGROUP_READEXON_VALUE_BOOL,iRegroupValue=REGROUP_READEXON_VALUE,fRegroupRatio=REGROUP_READEXON_THRESHOLD):
        
        #global REGROUP_READEXON_VALUE_BOOL
        #REGROUP_READEXON_VALUE_BOOL=bRegroupByValue
        #global REGROUP_READEXON_VALUE
        #REGROUP_READEXON_VALUE=iRegroupValue
        #global REGROUP_READEXON_THRESHOLD
        #REGROUP_READEXON_THRESHOLD=fRegroupRatio
        
        self.gene_start=None
        self.gene_end=None
        self.gene_strand=None
        self.size=None
        self.readGap=None
        self.transcript_to_exon={}
        self.set_emptyExonList()
        self.group_exon=[]
        
        self.parse_pafFile(sRefFile,sGeneId)
        
        self.sort_exonList()
        self.assign_exonId()
        self.update_transcriptList()
        
        #self.describe_transcript()
        #Temp=self.describe_exon()
        #print(Temp)
        #print("============================")
        self.fusion_firstReadStart()
        self.fusion_LastReadStop()
        if len(self.get_exonList())>1:
            self.fusion_CommonReadExon()
        self.sort_exonList()
        self.assign_exonId()
        #self.describe_transcript()
        #Temp=self.describe_exon()
        #print(Temp)
        
        self.make_group_exon()
        #sTemp=self.describe_group_exon(True)
        #print(sTemp)
        
    def compare_group_exon_with(self,oTranscriptContent):
        tAllGroup=self.get_groupExon()
        for iIndex in range(len(tAllGroup)):
            tGroup=tAllGroup[iIndex]
            for oSelfExon in tGroup:
                iMinVal=None
                iMaxVal=None
                for oExon in tGroup:
                    if iMinVal is None:
                        iMinVal=oExon.get_start()
                        iMaxVal=oExon.get_stop()
                    else:
                        iMinVal=min(iMinVal,oExon.get_start())
                        iMaxVal=max(iMaxVal,oExon.get_stop())
            tCorrespondingTrExon=[]
            for oTrExon in oTranscriptContent.get_exonList():
                if oTrExon.get_start()>iMaxVal or oTrExon.get_stop()<iMinVal:
                    #No overlap, pass
                    continue
                tCorrespondingTrExon.append(oTrExon)
            print("Group {} ({}-{}, {} elements) corresponding to :".format(iIndex,iMinVal,iMaxVal,len(tGroup)))
            if len(tCorrespondingTrExon)>0:
                for oTrExon in tCorrespondingTrExon:
                    print(oTrExon.get_index(),oTrExon.get_start(),oTrExon.get_stop())
            else:
                print("None")
                
    
    def describe_group_exon(self,resume=False):
        bResume=resume
        sContent=""
        tAllGroup=self.get_groupExon()
        if bResume:
            sContent+="GroupId\tVariant\tStart\tEnd\n"
        for iIndex in range(len(tAllGroup)):
            tGroup=tAllGroup[iIndex]
            if not bResume:
                sContent+="Group {} :\n".format(iIndex)
                for oExon in tGroup:
                    sContent+="{}\t{}\t{}\n".format(oExon.get_index(),oExon.get_start(),oExon.get_stop())
            else:
                sContent+="{}\t{}".format(iIndex,len(tGroup))
                iMinVal=None
                iMaxVal=None
                for oExon in tGroup:
                    if iMinVal is None:
                        iMinVal=oExon.get_start()
                        iMaxVal=oExon.get_stop()
                    else:
                        iMinVal=min(iMinVal,oExon.get_start())
                        iMaxVal=max(iMaxVal,oExon.get_stop())
                sContent+="\t{}\t{}\n".format(iMinVal,iMaxVal)
        return sContent
    
    def make_group_exon(self):
        self.set_emptyGroupExon()
        tGroupOfExon=[]
        tListOfExon=sorted(self.get_exonList(), key=lambda x: x.start)
        #tMadeExon=[]
        #tTodoExon=[]
        while len(tListOfExon)!=0:
            #if len(tTodoExon)==0:
                #tTodoExon.append(tListOfExon.pop())
            tMadeExon=[]
            tTodoExon=[tListOfExon.pop(0)]
            while len(tTodoExon)!=0:
                oTargetedExon=tTodoExon.pop(0)
                for oExon in tListOfExon:
                    if oExon.get_start()>oTargetedExon.get_stop() or oExon.get_stop()<oTargetedExon.get_start():
                        #No overlap, pass
                        continue
                    tTodoExon.append(oExon)
                for oExon in tTodoExon:
                    if oExon in tListOfExon:
                        tListOfExon.remove(oExon)
                tMadeExon.append(oTargetedExon)
            tOneGroupOfExon=[X for X in sorted(tMadeExon, key=lambda x: x.start)]
            tGroupOfExon.append(tOneGroupOfExon)        
        self.set_groupExon(tGroupOfExon)
        

    def set_emptyGroupExon(self):
        self.group_exon=[]
    
    def get_groupExon(self):
        return self.group_exon
        
    def set_groupExon(self,tList):
        self.group_exon=list(tList)
    
    def fusion_CommonReadExon(self):
        self.sort_exonList()
        
        #Regroup similar Exon
        dExon2SimilarExon={}
        tListOfExon=self.get_exonList()
        for iCurrentIndex in range(len(tListOfExon)):
            oCurrentExon=tListOfExon[iCurrentIndex]
            dExon2SimilarExon[oCurrentExon]=[oCurrentExon]
            iCurrentStart=oCurrentExon.get_start()
            iCurrentStop=oCurrentExon.get_stop()
            iCurrentSize=iCurrentStop-iCurrentStart+1
            for iAnotherIndex in range(len(tListOfExon)):
                if iAnotherIndex==iCurrentIndex:
                    continue
                oAnotherExon=tListOfExon[iAnotherIndex]
                iAnotherStart=oAnotherExon.get_start()
                iAnotherStop=oAnotherExon.get_stop()
                iAnotherSize=iAnotherStop-iAnotherStart+1
                iRefSize=max(iAnotherSize,iCurrentSize)
                if REGROUP_READEXON_VALUE_BOOL:
                    fThresholdValue=REGROUP_READEXON_VALUE
                else:
                    fThresholdValue=iRefSize*REGROUP_READEXON_THRESHOLD
                if iAnotherStop<iCurrentStart or iAnotherStart>iCurrentStop:
                    #Nothing in common, continue
                    continue
                elif iAnotherStop<=iCurrentStop and iAnotherStart>=iCurrentStart:
                    #Another is include into Current
                    if iCurrentStart+fThresholdValue>=iAnotherStart \
                    and iCurrentStop-fThresholdValue<=iAnotherStop:
                        #low distance between same extremity, they are maybe common
                        dExon2SimilarExon[oCurrentExon].append(oAnotherExon)
                    else:
                        #large distance between same extremity
                        continue
                elif iAnotherStop>=iCurrentStop and iAnotherStart<=iCurrentStart:
                    #Current is include into Another
                    if iCurrentStart-fThresholdValue<=iAnotherStart \
                    and iCurrentStop+fThresholdValue>=iAnotherStop:
                        #low distance between same extremity, they are maybe common
                        dExon2SimilarExon[oCurrentExon].append(oAnotherExon)
                    else:
                        #large distance between same extremity
                        continue
                elif iAnotherStop<=iCurrentStop and iAnotherStart<=iCurrentStart:
                    #Another overlap the beginning of Current
                    if iCurrentStart-fThresholdValue<=iAnotherStart \
                    and iCurrentStop-fThresholdValue<=iAnotherStop:
                        #low distance between same extremity, they are maybe common
                        dExon2SimilarExon[oCurrentExon].append(oAnotherExon)
                    else:
                        #large distance between same extremity
                        continue
                elif iAnotherStart>iCurrentStart and iAnotherStop>iCurrentStop:
                    #Another overlap the ending of Current
                    if iCurrentStop+fThresholdValue>=iAnotherStop \
                    and iCurrentStart+fThresholdValue>=iAnotherStart:
                        #low distance between same extremity, they are maybe common
                        dExon2SimilarExon[oCurrentExon].append(oAnotherExon)
                    else:
                        #large distance between same extremity
                        continue

        #for oCurrentExon in dExon2SimilarExon:
            #print(oCurrentExon.get_index(),":")
            #print([X.get_index() for X in dExon2SimilarExon[oCurrentExon]])
        #print("---------------")
        
        #Merge common exon
        tUsedExon=[]
        for oCurrentExon in dExon2SimilarExon:
            #print("===================")
            #print(oCurrentExon.get_index())
            #print("===================")
            if oCurrentExon not in tUsedExon:
                #Get all common
                tUsedExon.append(oCurrentExon)
                tSimilarExon=[oCurrentExon]
                tSimilarExon.extend([X for X in dExon2SimilarExon[tSimilarExon[0]] if X not in tUsedExon])
                #print(">>",tSimilarExon[0].get_index())
                #print("->",[X.get_index() for X in dExon2SimilarExon[tSimilarExon[0]] if X not in tUsedExon])
                bNoChange=False
                while not bNoChange:
                    iListSize=len(tSimilarExon)
                    for iIndex in range(iListSize):
                        #print(">>",tSimilarExon[iIndex].get_index())
                        #print("->",[X.get_index() for X in dExon2SimilarExon[tSimilarExon[iIndex]] if X not in tUsedExon])
                        tTargetedExon=[X for X in dExon2SimilarExon[tSimilarExon[iIndex]] if X not in tUsedExon]
                        tSimilarExon.extend(tTargetedExon)
                        tUsedExon.extend(tTargetedExon)
                    tSimilarExon=list(set(tSimilarExon))
                    tUsedExon=list(set(tUsedExon))
                    if len(tSimilarExon)==iListSize:
                        bNoChange=True
                #print("*****************Used*****************")
                #print([X.get_index() for X in tUsedExon])
                #print("*****************Similar*****************")
                #print([X.get_index() for X in tSimilarExon])
                #print("*****************")
                #Assign Start/Stop by majority
                dStop={}
                dStart={}
                for oSimilarExon in tSimilarExon:
                    try:
                        dStart[oSimilarExon.get_start()]+=1
                    except KeyError:
                        dStart[oSimilarExon.get_start()]=1
                    try:
                        dStop[oSimilarExon.get_stop()]+=1
                    except KeyError:
                        dStop[oSimilarExon.get_stop()]=1
                iMaxStart=None
                bStartEquality=False
                for iStartValue in dStart:
                    if iMaxStart is None:
                        iMaxStart=iStartValue
                    else:
                        if dStart[iMaxStart]==dStart[iStartValue]:
                            bEquality=True
                        elif dStart[iMaxStart]<dStart[iStartValue]:
                            iMaxStart=iStartValue
                if bStartEquality:
                    exit("ERROR 1807 : Start equality")
                iMaxStop=None
                bStopEquality=False
                for iStopValue in dStop:
                    if iMaxStop is None:
                        iMaxStop=iStopValue
                    else:
                        if dStop[iMaxStop]==dStop[iStopValue]:
                            bEquality=True
                        elif dStop[iMaxStop]<dStop[iStopValue]:
                            iMaxStop=iStopValue
                if bStopEquality:
                    exit("ERROR 1819 : Stop equality")
                oNewExon=ExonContent(iMaxStart,iMaxStop,"ThisIsNOTaTranscriptId")
                oNewExon.set_emptyTranscript()
                for oSimilarExon in tSimilarExon:
                    #print(oSimilarExon.describe_self())
                    oNewExon.update(oSimilarExon)
                    self.remove_exon(oSimilarExon)
                self.add_exon(oNewExon)
                
                
        
    
    def fusion_firstReadStart(self):
        self.sort_exonList()
        tTargetedExon=[]
        tTargetedExonIndex=[]
        tListOfTranscript=list(self.get_transcript())
        
        #get all first exon index of transcript
        for sTranscriptId in tListOfTranscript:
            tTargetedExonIndex.append(self.get_transcriptRelation(sTranscriptId)[0])
        tTargetedExonIndex=list(set(tTargetedExonIndex))

        #get all first exon of transcript, store quantity value for all stop
        dStop2Exon={}
        for oExon in self.get_exonList():
            if oExon.get_index() in tTargetedExonIndex:
                tTargetedExonIndex.remove(oExon.get_index())
                tTargetedExon.append(oExon)
                try:
                    dStop2Exon[oExon.get_stop()].append(oExon)
                except KeyError:
                    dStop2Exon[oExon.get_stop()]=[oExon]
            if len(tTargetedExonIndex)==0:
                break
                
        #Reorder tTargetedExon by stop population
        dQuantity2Stop={}
        tReorderedTargetedExon=[]
        for iStopValue in dStop2Exon:
            try:
                dQuantity2Stop[len(dStop2Exon[iStopValue])].append(iStopValue)
            except KeyError:
                dQuantity2Stop[len(dStop2Exon[iStopValue])]=[iStopValue]
        for iQuantity in sorted(dQuantity2Stop.keys(),reverse=True):
            for iTargetedStop in dQuantity2Stop[iQuantity]:
                for oExon in sorted(dStop2Exon[iTargetedStop], key=lambda x: x.start):
                    tReorderedTargetedExon.append(oExon)
        tTargetedExon=list(tReorderedTargetedExon)

        #fusion if they have the same stop +/-10%
        fThreshold=REGROUP_FIRSTREADEXON_THRESHOLD
        bNoChange=False
        while not bNoChange:
            bNoChange=True
            for iCurrentIndex in range(len(tTargetedExon)-1):
                oRefExon=tTargetedExon[iCurrentIndex]
                iRefExonSize=oRefExon.get_stop()-oRefExon.get_start()+1
                fFraction=iRefExonSize*fThreshold
                for iAnotherIndex in range(iCurrentIndex+1,len(tTargetedExon)):
                    oAnotherExon=tTargetedExon[iAnotherIndex]
                    if oAnotherExon.get_stop()>=oRefExon.get_stop()-fFraction and oAnotherExon.get_stop()<=oRefExon.get_stop()+fFraction:
                        oRefExon.update(oAnotherExon)
                        self.remove_exon(oAnotherExon)
                        bNoChange=False
                        break
                if not bNoChange:
                    tTargetedExon.remove(oAnotherExon)
                    break

    def fusion_LastReadStop(self):
        self.sort_exonList()
        tTargetedExon=[]
        tTargetedExonIndex=[]
        tListOfTranscript=list(self.get_transcript())
        
        #get all last exon index of transcript
        for sTranscriptId in tListOfTranscript:
            tTargetedExonIndex.append(self.get_transcriptRelation(sTranscriptId)[-1])
        tTargetedExonIndex=list(set(tTargetedExonIndex))
        
        #get all last exon of transcript
        dStart2Exon={}
        for oExon in self.get_exonList()[::-1]:
            if oExon.get_index() in tTargetedExonIndex:
                tTargetedExonIndex.remove(oExon.get_index())
                tTargetedExon.append(oExon)
                try:
                    dStart2Exon[oExon.get_start()].append(oExon)
                except KeyError:
                    dStart2Exon[oExon.get_start()]=[oExon]
            if len(tTargetedExonIndex)==0:
                break
                
        #Reorder tTargetedExon by start population
        dQuantity2Start={}
        tReorderedTargetedExon=[]
        for iStopValue in dStart2Exon:
            try:
                dQuantity2Start[len(dStart2Exon[iStopValue])].append(iStopValue)
            except KeyError:
                dQuantity2Start[len(dStart2Exon[iStopValue])]=[iStopValue]
        for iQuantity in sorted(dQuantity2Start.keys(),reverse=True):
            for iTargetedStart in dQuantity2Start[iQuantity]:
                for oExon in sorted(dStart2Exon[iTargetedStart],reverse=True, key=lambda x: x.stop):
                    tReorderedTargetedExon.append(oExon)
        tTargetedExon=list(tReorderedTargetedExon)        

        #fusion if they have the same start +/-10%
        fThreshold=REGROUP_LASTREADEXON_THRESHOLD
        bNoChange=False
        while not bNoChange:
            bNoChange=True
            for iCurrentIndex in range(len(tTargetedExon)-1):
                oRefExon=tTargetedExon[iCurrentIndex]
                iRefExonSize=oRefExon.get_stop()-oRefExon.get_start()+1
                fFraction=iRefExonSize*fThreshold
                for iAnotherIndex in range(iCurrentIndex+1,len(tTargetedExon)):
                    oAnotherExon=tTargetedExon[iAnotherIndex]
                    if oAnotherExon.get_start()>=oRefExon.get_start()-fFraction and oAnotherExon.get_start()<=oRefExon.get_start()+fFraction:
                        oRefExon.update(oAnotherExon)
                        self.remove_exon(oAnotherExon)
                        bNoChange=False
                        break
                if not bNoChange:
                    tTargetedExon.remove(oAnotherExon)
                    break
                    

    def parse_pafFile(self,sPathFile,sGeneId):
        
        dRead2Coord={}
        for sLine in open(sPathFile):
            sLine=sLine.strip()
            tLine=sLine.split()
            #tRef=["ReadId","ReadSize","ReadStart","ReadEnd","Strand","GeneId","GeneSize","GeneStart","GeneStop","NbrNotIndel","AlignSize","Quality"]
            dRef2Content={}
            for iIndex in range(len(PAF_COLUMN)):
                try:
                    dRef2Content[PAF_COLUMN[iIndex]]=int(tLine[iIndex])
                except ValueError:
                    dRef2Content[PAF_COLUMN[iIndex]]=tLine[iIndex]
            sCurrentTrId=dRef2Content["ReadId"]
            if sCurrentTrId not in self.get_transcript():
                self.add_transcript(sCurrentTrId)

            try:
                dRead2Coord[sCurrentTrId][min(dRef2Content["ReadStart"],dRef2Content["ReadStop"])]=(
                    dRef2Content["ReadStart"],
                    dRef2Content["ReadStop"],
                    dRef2Content["GeneStart"],
                    dRef2Content["GeneStop"],
                    #min(dRef2Content["ReadStart"],dRef2Content["ReadStop"]),
                    #max(dRef2Content["ReadStart"],dRef2Content["ReadStop"]),
                    #min(dRef2Content["GeneStart"],dRef2Content["GeneStop"]),
                    #max(dRef2Content["GeneStart"],dRef2Content["GeneStop"]),
                    dRef2Content["Strand"]
                    )
            except KeyError:
                dRead2Coord[sCurrentTrId]={
                    min(dRef2Content["ReadStart"],dRef2Content["ReadStop"]):
                        (
                            dRef2Content["ReadStart"],
                            dRef2Content["ReadStop"],
                            dRef2Content["GeneStart"],
                            dRef2Content["GeneStop"],
                            #min(dRef2Content["ReadStart"],dRef2Content["ReadStop"]),
                            #max(dRef2Content["ReadStart"],dRef2Content["ReadStop"]),
                            #min(dRef2Content["GeneStart"],dRef2Content["GeneStop"]),
                            #max(dRef2Content["GeneStart"],dRef2Content["GeneStop"]),
                            dRef2Content["Strand"]
                        )
                    }
            oExon=ExonContent(int(dRef2Content["GeneStart"]),int(dRef2Content["GeneStop"]),sCurrentTrId)
            self.update_exonList(oExon)
        self.set_geneCoord(int(dRef2Content["GeneStart"]),int(dRef2Content["GeneStop"]),1)
        self.set_size(int(dRef2Content["ReadSize"]))
        self.compute_readGap(dRead2Coord)
        
    def get_size(self):
        return self.size
        
    def set_size(self,iInt):
        self.size=iInt

    def compute_readGap(self,dDict):
        dRead2GapCoord={}
        for sTrId in dDict:
            iLastValue=None
            iLastLimit=None
            #print(dDict[sTrId])
            for iCoordTag in sorted(dDict[sTrId]):
                iReadStart=dDict[sTrId][iCoordTag][0]
                iReadStop=dDict[sTrId][iCoordTag][1]
                iGeneStart=dDict[sTrId][iCoordTag][2]
                iGeneStop=dDict[sTrId][iCoordTag][3]
                sReadStrand=dDict[sTrId][iCoordTag][4]
                #print("===================")
                #print(iReadStart,iReadStop)
                #print(iGeneStart,iGeneStop)
                #print(sReadStrand)
                if iLastValue==None:
                    iLastValue=iReadStop
                    if sReadStrand=="+":
                        iLastLimit=iGeneStop
                    else:
                        iLastLimit=iGeneStart
                else:
                    iSpaceBetween=(iReadStart-1)-iLastValue
                    
                    if iSpaceBetween>0:
                        #There is a gap
                        dbRegion=(iLastLimit,iGeneStop)
                        #print(iSpaceBetween)
                        #print(dbRegion)
                        #input()
                        if sReadStrand=="+":
                            try:
                                dRead2GapCoord[sTrId][iLastLimit]=(iLastLimit,iGeneStart,iSpaceBetween)
                            except KeyError:
                                dRead2GapCoord[sTrId]={iLastLimit:(iLastLimit,iGeneStart,iSpaceBetween)}
                        else:
                            try:
                                dRead2GapCoord[sTrId][iGeneStart]=(iGeneStart,iLastLimit,iSpaceBetween)
                            except KeyError:
                                dRead2GapCoord[sTrId]={iGeneStart:(iGeneStart,iLastLimit,iSpaceBetween)}
                    if iLastValue<iReadStop:
                        iLastValue=iReadStop
                        if sReadStrand=="+":
                            iLastLimit=iGeneStop
                        else:
                            iLastLimit=iGeneStart
        self.set_readGap(dRead2GapCoord)
        
    def set_readGap(self,dDict):
        self.readGap=dDict
        #print(dDict)
        #exit()
        
    def get_readGap(self):
        return self.readGap
        
        
class ExonContent:

    def __init__(self,iCoordStart,iCoordStop,sTranscriptId):
        self.exon_id=None
        self.start=iCoordStart
        self.stop=iCoordStop
        self.size=iCoordStop-iCoordStart+1
        self.transcript=[sTranscriptId]
    
    def get_index(self):
        return self.exon_id
    
    def get_start(self):
        return self.start
    
    def set_start(self,iValue):
        self.start=iValue
    
    def get_stop(self):
        return self.stop
        
    def set_stop(self,iValue):
        self.stop=iValue
        
    def get_transcript(self):
        return self.transcript
    
    def update(self,oObject):
        self.add_transcript(oObject.get_transcript())
        
    def add_transcript(self,tList):
        for oObject in tList:
            self.transcript.append(oObject)

    def set_emptyTranscript(self):
        self.transcript=[]

    def set_index(self,oValue):
        self.exon_id=oValue

    def describe_self(self,bStdout=False):
        sContent=""
        for oValue in [self.exon_id,self.start,self.stop,self.size,self.transcript]:
            if sContent!="":
                sContent+="\t"
            if isinstance(oValue, (list, tuple)):
                sContent+=",".join(oValue)
            else:
                sContent+="{}".format(oValue)
        if bStdout:
            print(sContent)
        sContent+="\n"
        return sContent
