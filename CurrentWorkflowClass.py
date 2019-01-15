# coding: utf-8
"""Python3.6"""

import os
import sys
import re
import copy
import random
import subprocess

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

CORRECT_ALIGN_STEP__EXTAND_VALUE=0
LALIGN_LAUNCHER="~/Git/fasta36/bin/lalign36 -3 -U" #WARNING: only forwad strand
#LALIGN_THRESHOLD=80

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
def GetReverseComplement(sString):
    #print(sString)
    #print("".join([REVERSE_COMP[X] for X in sString[::-1]]))
    return "".join([REVERSE_COMP[X] for X in sString[::-1]])

def ParseLalignResult(sFile,sReadBlock):
    bTrackingAlignment=False
    iAlignmentIndex=0
    sReadString=""
    sGeneString=""
    sAlignString=""
    bStop=False
    for sLine in open(sFile):
        sLine=sLine.strip()
        if ">>>" in sLine:
            tLine=sLine.split(" ")
            iReadSize=int(tLine[-3])
        elif ">>" in sLine and not bStop:
            tLine=re.split("[>< ()]",sLine.replace(">>",""))
            sGeneBlock=tLine[0]
            iGeneSize=int(tLine[-3])
            bStop=True
        elif "!! No sequences" in sLine:
            return {}
            #exit("No sequence")
        elif "identity" in sLine:
            tLine=re.split("[%>< ():-]",sLine)
            iReadStart=int(tLine[-5])
            iReadStop=int(tLine[-4])
            iGeneStart=int(tLine[-3])
            iGeneStop=int(tLine[-2])
            iAlignSize=int(tLine[-9])
            bTrackingAlignment=True
            continue
        elif "query   sequences" in sLine or (">>" in sLine and bStop) or (">--" in sLine and bStop):
            bTrackingAlignment=False
            try:
                oCrashMe=iReadStart
            except UnboundLocalError:
                return {}
            if iReadStart>iReadStop:
                sStrand="-"
                iReadStart, iReadStop = iReadStop, iReadStart
            else:
                sStrand="+"
            iAlignGap=sReadString.count("-")+sGeneString.count("-")
            iAlignIdentity=iAlignSize-sAlignString.count(" ")
            return {
                    "ReadId":sReadBlock,
                    "ReadSize":iReadSize,
                    "ReadStart":iReadStart,
                    "ReadStop":iReadStop,
                    "Strand":sStrand,
                    "GeneId":sGeneBlock,
                    "GeneSize":iGeneSize,
                    "GeneStart":iGeneStart,
                    "GeneStop":iGeneStop,
                    #"AlignIdentity":iAlignIdentity,
                    "AlignIdentity":len([X for X in range(len(sReadString)) if sReadString[X]==sGeneString[X]])/len(sReadString),
                    "AlignGap":iAlignGap,
                    "AlignSize":iAlignSize,
                    "GeneString":sGeneString,
                    "ReadString":sReadString,
                    "AlignString":sAlignString
                    }
        elif bTrackingAlignment:
            tLine=sLine.split(" ")
            if iAlignmentIndex==0:
                if tLine[0]==sReadBlock:
                    sReadString+=tLine[-1]
                    iAlignmentIndex=1
            elif iAlignmentIndex==1:
                #print("____",sLine,"_____")
                sAlignString+=sLine.replace(":","|")
                iAlignmentIndex=2
            elif iAlignmentIndex==2:
                sGeneString+=tLine[-1]
                iAlignmentIndex=0

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

def WriteTempFasta(dDict):
    sName="TempFasta"+str(random.random())+".fa"
    FileTemp=open(sName,"w")
    for sKey in dDict:
        FileTemp.write(">"+sKey+"\n"+dDict[sKey]+"\n")
    FileTemp.close()
    return sName

def ExecuteBashCommand(scriptfile,bFirstLine=True):
    #print("DEBUG:{}".format(scriptfile))
    #return None
    sName="TempCommand"+str(random.random())+".sh"
    FileTemp=open(sName,"w")
    FileTemp.write(scriptfile)
    FileTemp.close()

    tTable=LaunchBashFile(["bash",sName])

    if tTable[0]==0:
        if tTable[1]:
            print(tTable[1].decode("utf-8"))
        print("Bash command done")
    else:
        print(tTable[2])
        print("Bash command fail: {}\nScript bash : {}".format(scriptfile,sName))
    
    os.remove(sName)

def LaunchBashFile(cmdArray):
	sp = subprocess.Popen(cmdArray, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = sp.communicate()
	if err:
	    print("standard error of subprocess:")
	    print(err.decode("utf-8"))
	    #if out:
            #print("standard output of subprocess:")
            #print(out.decode("utf-8"))
	print("returncode of subprocess:", sp.returncode)
	return [sp.returncode,out,err]

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
    
    def __init__(self,sGeneId,sReadXml,sOutputFile,sFastaFile,sRefFile):
        
        #sReadPaf=sReadXml.replace("xml","paf")
        
        self.size=None
        self.gene_id=sGeneId
        
        self.tr_dict_seq=Fasta2Dict(sFastaFile)
        self.ref_seq=sorted(Fasta2Dict(sRefFile).values())[0]
        
        self.currentMatrix=[]
        self.matrix_lineName=[]
        self.globalVector=[]
        self.globalBlockName=[]
        self.currentPopVector=[]
        self.globalPopVector=[]
        self.blockConfidence={}
        
        self.tr_dict_data={}
        self.tr_dict_data_alt={}
        self.parse_xmlFile(sReadXml,sGeneId)
        #print(self.get_tr_dict_data())
        #self.parse_pafFile(sReadPaf,sGeneId)
        #print(self.get_tr_dict_data())
        ##exit()
        self.setup_matrix()
        #print("DebugMe:",self.get_submatrix_byLine(2,2)[0][13357:13439])
        
        #print(self.get_submatrix_byLine(2,2))
        #exit()
        
        self.setup_popVector()
        self.setup_globalVector()
        self.assign_blockName()
        self.assign_blockConfidence()
        
        #self.globalData2tsv("PrimaryTable.tsv")
        
        self.regroup_globalVector()  #This will modify the matrix!!
        #print("DebugMe:",self.get_submatrix_byLine(2,2)[0][13357:13439])
        
        self.setup_popVector()
        self.setup_globalVector()
        self.assign_blockName()
        self.assign_blockConfidence()

        #print("selfBlock:",self.get_limitedBlockName())
        #print("selfBlockSize:",self.get_blockSize())
        #print("selfBlockCoord:",self.get_blockCoord())
        #print("selfBlockPop:",self.get_blockPop())
            
        self.globalData2tsv(sOutputFile)
    
    #def debug_matrix(self):
        #dTr2Data=self.get_tr_dict_data()
        #tTrId=self.get_matrix_lineName()
        #for iTrId in range(len(tTrId)):
            #sTrId=tTrId[iTrId]
            #print("~~~~~~~~~~~~~~~~~")
            #print(sTrId,dTr2Data[sTrId][sorted(dTr2Data[sTrId])[0]]["Strand"])
            
            #tCurrentVector=self.get_submatrix_byLine(iTrId,iTrId)[0]
                        
            #for iDbIndex in range(len(dTr2Data[sTrId])):
                #dbCurrentKey=sorted(dTr2Data[sTrId])[iDbIndex]
                #iCurrentGeneStart=dTr2Data[sTrId][dbCurrentKey]["GeneStart"]
                #iCurrentGeneStop=dTr2Data[sTrId][dbCurrentKey]["GeneStop"]
                
                #print("------------")
                #print("TrData",dTr2Data[sTrId][dbCurrentKey])
                #print(tCurrentVector[iCurrentGeneStart:iCurrentGeneStop+1])
                
                #if sum(tCurrentVector[iCurrentGeneStart:iCurrentGeneStop+1])!=len(tCurrentVector[iCurrentGeneStart:iCurrentGeneStop+1]):
                    #exit("Error 1872 : Vector not correcting filled")
        

    
    def get_geneId(self):
        return self.gene_id
    
    def get_ref_seq(self,dbCoord=False):
        if dbCoord:
            return self.ref_seq[dbCoord[0]:dbCoord[1]+1]
        else:
            return self.ref_seq
    
    def get_read_seq(self,sReadName,dbCoord=False):
        if dbCoord:
            return self.tr_dict_seq[sReadName][dbCoord[0]:dbCoord[1]:+1]
        else:
            return self.tr_dict_seq[sReadName]
        
    def print_read_seq(self):
        for iLineIndex in range(len(self.get_matrix())):
            tLineModel=self.get_line_BlockName_Structure(iLineIndex)
            print(">"+self.get_matrix_lineName(iLineIndex)+"\n"+self.get_read_seq(self.get_matrix_lineName(iLineIndex)))
    
    def apply_correction(self):
        #print("DebugMe:",self.get_submatrix_byLine(2,2)[0][13357:13439])
        bMatrixIsModified=True
        iCorrectionStep=0
        while bMatrixIsModified:
            iCorrectionStep+=1
            bMatrixIsModified=False
            if self.correct_misalignment_onref():
                bMatrixIsModified=True

            if bMatrixIsModified:
                self.regroup_globalVector()
                self.setup_popVector()
                self.setup_globalVector()
                self.assign_blockName()
                self.assign_blockConfidence()
                
                self.globalData2tsv("Intermediate{}.spliceSummary.tsv".format(iCorrectionStep))
                #print("DebugMe:",self.get_submatrix_byLine(2,2)[0][13357:13439])
    
    def make_groupOfBlock(self,tListOfBlock):
        tGroupOfGroupOfBlock=[]
        tCurrentGroupOfBlock=[]
        for sElement in tListOfBlock:
            if sElement=="i":
                if len(tCurrentGroupOfBlock)!=0:
                    tGroupOfGroupOfBlock.append(list(tCurrentGroupOfBlock))
                tCurrentGroupOfBlock=[]
            elif sElement!=tListOfBlock[-1]:
                tCurrentGroupOfBlock.append(sElement)
            else:
                tCurrentGroupOfBlock.append(sElement)
                tGroupOfGroupOfBlock.append(list(tCurrentGroupOfBlock))
        return tGroupOfGroupOfBlock
    
    def get_suspiciousblock(self,tGroupOfBlock,tGroupOfGroupOfBlock):
        tSuspiciousBlock=[]
        for iIndex in range(len(tGroupOfGroupOfBlock)):
            tCurrentGroupOfBlock=tGroupOfGroupOfBlock[iIndex]
            tCurrentPop=[self.get_blockPop(X) for X in tCurrentGroupOfBlock]
            iMaxPop=max(tCurrentPop)
                        
            if iIndex==0:
                if iMaxPop==tCurrentPop[-1]:
                    continue
                else:
                    tCurrentPop=tCurrentPop[-(tCurrentPop[::-1].index(iMaxPop)+1):]
                    tCurrentGroupOfBlock[:-len(tCurrentPop)]
            elif iIndex==len(tGroupOfGroupOfBlock)-1:
                if iMaxPop==tCurrentPop[0]:
                    continue
                else:
                    tCurrentPop=tCurrentPop[:tCurrentPop.index(iMaxPop)+1]
                    tCurrentGroupOfBlock[-len(tCurrentPop):]
            tSuspiciousBlock+=[tCurrentGroupOfBlock[X] for X in range(len(tCurrentPop)) if tCurrentPop[X]==1 or tCurrentPop[X]<=0.1*iMaxPop]
        return tSuspiciousBlock
    
    def get_PreviousConfidentBlock(self,tSuspiciousBlock,tCurrentSuspiciousBlock,sCurrentBlock,tReadBlock,tGroupOfReadBlock):
        sPreviousBlock=None
        if sCurrentBlock==tGroupOfReadBlock[0][0]:
            sPreviousBlock=None
        if len(tGroupOfReadBlock[0])>1:
            if sCurrentBlock==tGroupOfReadBlock[0][1] and isinstance(tGroupOfReadBlock[0][0],int):
                sPreviousBlock=None
        iPrevious=-1
        while sPreviousBlock is None:
            if tReadBlock.index(sCurrentBlock)+iPrevious<0:
                break
            sCandidat=tReadBlock[tReadBlock.index(sCurrentBlock)+iPrevious]
            #print("sCandidat",sCandidat)
            if sCandidat!="i" and sCandidat not in tSuspiciousBlock and not isinstance(sCandidat,int):
                sPreviousBlock=sCandidat
                break
            iPrevious-=1
        return sPreviousBlock
    
    def get_NextConfidentBlock(self,tSuspiciousBlock,tCurrentSuspiciousBlock,sCurrentBlock,tReadBlock,tGroupOfReadBlock):
        sNextBlock=None
        if sCurrentBlock==tGroupOfReadBlock[-1][-1]:
            sNextBlock=None
        if len(tGroupOfReadBlock[-1])>1:
            if sCurrentBlock==tGroupOfReadBlock[-1][-2] and isinstance(tGroupOfReadBlock[-1][-1],int):
                sNextBlock=None
        iNext=1
        while sNextBlock is None:
            try:
                sCandidat=tReadBlock[tReadBlock.index(sCurrentBlock)+iNext]
            except IndexError:
                break
            if sCandidat!="i" and sCandidat not in tSuspiciousBlock and not isinstance(sCandidat,int):
                sNextBlock=sCandidat
                break
            iNext+=1
        return sNextBlock
    
    def get_dbAlignGeneCoord(self,sLineName,sCurrentBlockCoord):
        tCorrespondingData=[X for X in self.get_tr_dict_data_alt()[sLineName].keys() if sCurrentBlockCoord[0] in range(X[0],X[1]+1) or sCurrentBlockCoord[1] in range(X[0],X[1]+1)]
        if len(tCorrespondingData)==0:
            print(self.get_tr_dict_data_alt()[sLineName].keys())
            exit("Error 1865 : no correspondance between Block and Read Alignment")
        elif len(tCorrespondingData)!=1:
            print(self.get_tr_dict_data_alt()[sLineName].keys())
            exit("Error 1868 : too much correspondance between Block and Read Alignment")
        return tCorrespondingData[0]
    
    #def get_alignReadSeq(self,sLineName,dbAlignGeneCoord,sCurrentBlockCoord,iExtand=0):
        #dTempDict=self.get_tr_dict_data_alt()[sLineName][dbAlignGeneCoord]
        
        #iGeneStart=dTempDict["GeneStart"]
        #iGeneStop=dTempDict["GeneStop"]
        #iReadStart=dTempDict["ReadStart"]
        #iReadStop=dTempDict["ReadStop"]
        #sGeneStrand=dTempDict["Strand"]
        #iDeltaStart=iGeneStart-sCurrentBlockCoord[0]
        #iDeltaStop=iGeneStop-sCurrentBlockCoord[1]
        #if sGeneStrand=="+":
            #iNewReadStart=iReadStart-iDeltaStart
            #iNewReadStop=iReadStop-iDeltaStop
        #else:
            #iNewReadStart=iReadStart+iDeltaStop
            #iNewReadStop=iReadStop+iDeltaStart
        #return self.get_read_seq(sLineName,(iNewReadStart-iExtand,iNewReadStop+iExtand))
    
    def extract_alignReadData(self,sLineName,dbAlignGeneCoord,dbCurrentBlockCoord):
        dTempDict=self.get_tr_dict_data_alt()[sLineName][dbAlignGeneCoord]
        sReadId=dTempDict["ReadId"]
        iReadSize=dTempDict["ReadSize"]
        sStrand=dTempDict["Strand"]
        sGeneId=dTempDict["GeneId"]
        iGeneSize=dTempDict["GeneSize"]
        sReadString=dTempDict["ReadString"]
        sGeneString=dTempDict["GeneString"]
        sAlignString=dTempDict["AlignString"]
        iGeneStart=dTempDict["GeneStart"]
        iGeneStop=dTempDict["GeneStop"]
        iReadStart=dTempDict["ReadStart"]
        iReadStop=dTempDict["ReadStop"]
        iDeltaStart=iGeneStart-dbCurrentBlockCoord[0]
        iDeltaStop=iGeneStop-dbCurrentBlockCoord[1]
        
        sUngappedReadString=sReadString.replace("-","")
        tUngappedReadString=list(range(iReadStart,iReadStop+1))
        tReadString=[]
        iCoord=iReadStart
        for iIndex in range(len(sReadString)):
            cChar=sReadString[iIndex]
            if cChar!="-":
                tReadString.append(iCoord)
                iCoord+=1
            else:
                tReadString.append(-1)
        
        sUngappedGeneString=sGeneString.replace("-","")
        tUngappedGeneString=list(range(iGeneStart,iGeneStop+1))
        tGeneString=[]
        iCoord=iGeneStart
        for iIndex in range(len(sGeneString)):
            cChar=sGeneString[iIndex]
            if cChar!="-":
                tGeneString.append(iCoord)
                iCoord+=1
            else:
                tGeneString.append(-1)
                
        #print(tGeneString)
        #print(dbCurrentBlockCoord)
        
        dData2Coord={}
        
        if dbCurrentBlockCoord[0]!=tGeneString[0]:
            iIndex=tGeneString.index(dbCurrentBlockCoord[0]-1)
            
            sUpGeneString=sGeneString[:iIndex+1]
            sUpReadString=sReadString[:iIndex+1]
            sUpAlignString=sAlignString[:iIndex+1]
            iUpReadStart=tReadString[0]
            iUpReadStop=tReadString[iIndex]
            iUpGeneStart=tGeneString[0]
            iUpGeneStop=tGeneString[iIndex]
            
            dData2Coord[(min(iUpGeneStart,iUpGeneStop),max(iUpGeneStart,iUpGeneStop))]={
                    "ReadId":sReadId,
                    "ReadSize":iReadSize,
                    "ReadStart":iUpReadStart,
                    "ReadStop":iUpReadStop,
                    "Strand":sStrand,
                    "GeneId":sGeneId,
                    "GeneSize":iGeneSize,
                    "GeneStart":iUpGeneStart,
                    "GeneStop":iUpGeneStop,
                    "AlignIdentity":len([X for X in range(len(sUpGeneString)) if sUpReadString[X]==sUpGeneString[X]])/len(sUpReadString),
                    "AlignGap":sUpGeneString.count("-")+sUpReadString.count("-"),
                    "AlignSize":len(sUpReadString),
                    "GeneString":sUpGeneString,
                    "ReadString":sUpReadString,
                    "AlignString":sUpAlignString
                    }
            print("Before",dData2Coord)
        
        if dbCurrentBlockCoord[-1]!=tGeneString[-1]:
            iIndex=tGeneString.index(dbCurrentBlockCoord[-1]+1)
            
            sDownGeneString=sGeneString[iIndex:]
            sDownReadString=sReadString[iIndex:]
            sDownAlignString=sAlignString[iIndex:]
            
            iDownReadStart=tReadString[iIndex]
            iDownReadStop=tReadString[-1]
            iDownGeneStart=tGeneString[iIndex]
            iDownGeneStop=tGeneString[-1]
            
            dData2Coord[(min(iDownGeneStart,iDownGeneStop),max(iDownGeneStart,iDownGeneStop))]={
                    "ReadId":sReadId,
                    "ReadSize":iReadSize,
                    "ReadStart":iDownReadStart,
                    "ReadStop":iDownReadStop,
                    "Strand":sStrand,
                    "GeneId":sGeneId,
                    "GeneSize":iGeneSize,
                    "GeneStart":iDownGeneStart,
                    "GeneStop":iDownGeneStop,
                    "AlignIdentity":len([X for X in range(len(sDownReadString)) if sDownReadString[X]==sDownGeneString[X]])/len(sDownReadString),
                    "AlignGap":sDownGeneString.count("-")+sDownReadString.count("-"),
                    "AlignSize":len(sDownReadString),
                    "GeneString":sDownGeneString,
                    "ReadString":sDownReadString,
                    "AlignString":sDownAlignString
                    }
            print("After",dData2Coord)
            
        iTargetIndexStart=tGeneString.index(dbCurrentBlockCoord[0])
        iTargetIndexStop=tGeneString.index(dbCurrentBlockCoord[-1])
        
        sTargetGeneString=sGeneString[iTargetIndexStart:iTargetIndexStop+1]
        sTargetReadString=sReadString[iTargetIndexStart:iTargetIndexStop+1]
        sTargetAlignString=sAlignString[iTargetIndexStart:iTargetIndexStop+1]
        iTargetReadStart=tReadString[iTargetIndexStart]
        iTargetReadStop=tReadString[iTargetIndexStop]
        iTargetGeneStart=tGeneString[iTargetIndexStart]
        iTargetGeneStop=tGeneString[iTargetIndexStop]
        
        dData2Coord[(min(iTargetGeneStart,iTargetGeneStop),max(iTargetGeneStart,iTargetGeneStop))]={
                "ReadId":sReadId,
                "ReadSize":iReadSize,
                "ReadStart":iTargetReadStart,
                "ReadStop":iTargetReadStop,
                "Strand":sStrand,
                "GeneId":sGeneId,
                "GeneSize":iGeneSize,
                "GeneStart":iTargetGeneStart,
                "GeneStop":iTargetGeneStop,
                "AlignIdentity":len([X for X in range(len(sTargetReadString)) if sTargetReadString[X]==sTargetGeneString[X]])/len(sTargetReadString),
                "AlignGap":sTargetGeneString.count("-")+sTargetReadString.count("-"),
                "AlignSize":len(sTargetReadString),
                "GeneString":sTargetGeneString,
                "ReadString":sTargetReadString,
                "AlignString":sTargetAlignString
                }
        
        print("All",dData2Coord)
        
        return (dData2Coord,(min(iTargetGeneStart,iTargetGeneStop),max(iTargetGeneStart,iTargetGeneStop)))
        
        
    def get_alignGeneSeq(self,tPotentialTarget,iExtand=0):
        dTarget2Seq={}
        for sTargetId in tPotentialTarget:
            dTarget2Seq[sTargetId]=self.get_ref_seq((self.get_blockCoord(sTargetId)[0]-iExtand,self.get_blockCoord(sTargetId)[1]+iExtand))
        return dTarget2Seq
    
    def correct_misalignment_onref(self):
        tGroupOfBlock=self.get_limitedBlockName()
        tGroupOfGroupOfBlock=self.make_groupOfBlock(tGroupOfBlock)
        print("tGroupOfGroupOfBlock",tGroupOfGroupOfBlock)
        
        tSuspiciousBlock=self.get_suspiciousblock(tGroupOfBlock,tGroupOfGroupOfBlock)
        if len(tSuspiciousBlock)==0:
            return False
        else:
            print("There is suspicious block : {}".format(tSuspiciousBlock))
        
        tGeneBlock=self.get_limitedBlockName()
        print("tGeneBlock",tGeneBlock)
        
        for iLineIndex in range(len(self.get_matrix())):
            sLineName=self.get_matrix_lineName(iLineIndex)
            tLineModel=self.get_line_BlockName_Structure(iLineIndex)
            tCurrentSuspiciousBlock=sorted(set(tLineModel) & set(tSuspiciousBlock))
                        
            if len(set(tLineModel) & set(tSuspiciousBlock))==0:
                continue
            else:
                print("{} have {} suspicious block on {}: {}".format(sLineName,len(tCurrentSuspiciousBlock),self.get_geneId(),tCurrentSuspiciousBlock))
            
            print(sLineName)
            tReadBlock=self.get_line_BlockName_Structure(iLineIndex)
            print("tReadBlock",tReadBlock)
            tGroupOfReadBlock=self.make_groupOfBlock(tReadBlock)
            print("tGroupOfReadBlock",tGroupOfReadBlock)
            for sCurrentBlock in tCurrentSuspiciousBlock:
                print(sCurrentBlock,self.get_blockCoord(sCurrentBlock))
                
                if sCurrentBlock==tReadBlock[0] or (sCurrentBlock==tReadBlock[1] and isinstance(tReadBlock[0],int)):
                    continue
                
                sPreviousBlock=self.get_PreviousConfidentBlock(tSuspiciousBlock,tCurrentSuspiciousBlock,sCurrentBlock,tReadBlock,tGroupOfReadBlock)
                sNextBlock=self.get_NextConfidentBlock(tSuspiciousBlock,tCurrentSuspiciousBlock,sCurrentBlock,tReadBlock,tGroupOfReadBlock)
                print("{0} chained into {1} {0} {2}".format(sCurrentBlock,sPreviousBlock,sNextBlock))
                print(self.get_alignGeneSeq([sCurrentBlock,sPreviousBlock,sNextBlock]))
                
                if sPreviousBlock is not None and sNextBlock is not None:
                    tPotentialTarget=[X for X in tGeneBlock[tGeneBlock.index(sPreviousBlock)+1:tGeneBlock.index(sNextBlock)] if X!="i" and not isinstance(X,int) and X not in tSuspiciousBlock]
                elif sPreviousBlock is None and sNextBlock is not None:
                    tPotentialTarget=[X for X in tGeneBlock[0:tGeneBlock.index(sNextBlock)] if X!="i" and not isinstance(X,int) and X not in tSuspiciousBlock]
                elif sPreviousBlock is not None and sNextBlock is None:
                    tPotentialTarget=[X for X in tGeneBlock[tGeneBlock.index(sPreviousBlock)+1:] if X!="i" and not isinstance(X,int) and X not in tSuspiciousBlock]
                else:
                    exit("Error 1904 : no confident block in the read")
                
                if len(tPotentialTarget)==0:
                    print("No potential target. Skip.")
                    continue
                
                print("Potential target are {}".format(tPotentialTarget))
                
                #dTarget2Seq=self.get_alignGeneSeq(tPotentialTarget,CORRECT_ALIGN_STEP__EXTAND_VALUE)
                dTarget2Seq=self.get_alignGeneSeq(tPotentialTarget)
                print(dTarget2Seq)
                
                sCurrentBlockCoord=self.get_blockCoord(sCurrentBlock)
                dbAlignGeneCoord=self.get_dbAlignGeneCoord(sLineName,sCurrentBlockCoord)
                #sAlignReadSeq=self.get_alignReadSeq(sLineName,dbAlignGeneCoord,sCurrentBlockCoord,CORRECT_ALIGN_STEP__EXTAND_VALUE)
                dbTemp=self.extract_alignReadData(sLineName,dbAlignGeneCoord,sCurrentBlockCoord)
                print(dbTemp)
                dDict=dbTemp[0]
                dbTargetSeq=dbTemp[1]
                sAlignReadSeq=dDict[dbTargetSeq]["ReadString"].replace("-","")
                print(sCurrentBlock,":",sAlignReadSeq)
                
                sTempFastaRead=WriteTempFasta({sCurrentBlock:sAlignReadSeq})
                sTempFastaRef=WriteTempFasta(dTarget2Seq)
                sTempResult="TempLALIGN"+str(random.random())+".out"
                
                ExecuteBashCommand("{0} {1} {2} > {3}".format(LALIGN_LAUNCHER,sTempFastaRead,sTempFastaRef,sTempResult))
                #dOldAlignData=self.get_tr_dict_data_alt()[sLineName][dbAlignGeneCoord]
                dOldAlignData=dDict[dbTargetSeq]
                
                #print(dOldAlignData)
                #print(dDict[dbTargetSeq])
                #exit()
                
                dNewAlignData=ParseLalignResult(sTempResult,sCurrentBlock)
                print("NewAlignData",dNewAlignData)
                
                #print(sTempFastaRead,sTempFastaRef,sTempResult)
                for sFile in [sTempFastaRead,sTempFastaRef,sTempResult]:
                    ExecuteBashCommand("rm {}".format(sFile))
                
                if len(dNewAlignData)==0:
                    print("No valuable alignment. Skip.")
                    continue
                
                fOldDistance=dOldAlignData["AlignIdentity"]
                fNewDistance=dNewAlignData["AlignIdentity"]
                
                if fNewDistance<fOldDistance:
                    exit("TODO2066: New alignment less good than old alignment")
                
                ##Check unaligned seq
                iUpstreamPenalty=dNewAlignData["ReadStart"]-1
                if iUpstreamPenalty<=REGROUP_GLOBALVECTOR_VALUE:
                    iUpstreamPenalty=None
                else:
                    #exit("ERROR 2239 : too long unaligned upstream to ignore it")
                    print("Warning 2239 : too long unaligned upstream to ignore it")
                iDownstreamPenalty=dNewAlignData["ReadSize"]-dNewAlignData["ReadStop"]
                if iDownstreamPenalty<=REGROUP_GLOBALVECTOR_VALUE:
                    iDownstreamPenalty=None
                else:
                    #exit("ERROR 2244 : too long unaligned downstream to ignore it")
                    print("Warning 2244 : too long unaligned downstream to ignore it")
                
                ##Update newAlignData with data from OldAlignData
                sAlignedGeneSeq=dTarget2Seq[dNewAlignData["GeneId"]]
                dbAlignedGeneCoord=self.get_blockCoord(dNewAlignData["GeneId"])
                
                dNewAlignData["ReadId"]=dOldAlignData["ReadId"]
                dNewAlignData["ReadSize"]=dOldAlignData["ReadSize"]
                dNewAlignData["GeneId"]=dOldAlignData["GeneId"]
                dNewAlignData["GeneSize"]=dOldAlignData["GeneSize"]
                dNewAlignData["Strand"]=dOldAlignData["Strand"]
                dNewAlignData["ReadStop"]=dDict[dbTargetSeq]["ReadStop"]-(len(sAlignReadSeq)-dNewAlignData["ReadStop"])
                dNewAlignData["ReadStart"]=dDict[dbTargetSeq]["ReadStart"]+dNewAlignData["ReadStart"]-1
                dNewAlignData["GeneStop"]=dbAlignedGeneCoord[1]-(len(sAlignedGeneSeq)-dNewAlignData["GeneStop"])
                dNewAlignData["GeneStart"]=dbAlignedGeneCoord[0]+dNewAlignData["GeneStart"]-1
                
                #DEBUG
                print("--------------Replace--------------")
                print("Old alignment distance : {} \n{}\n{}\n{}".format(fOldDistance,
                    dOldAlignData["GeneString"],dOldAlignData["AlignString"],dOldAlignData["ReadString"]))
                print("New alignment distance : {} \n{}\n{}\n{}".format(fNewDistance,
                    dNewAlignData["GeneString"],dNewAlignData["AlignString"],dNewAlignData["ReadString"]))
                print("--------------/Replace--------------")
                #/DEBUG
                
                print(dDict)
                print(dbTargetSeq)
                #exit()
                
                del dDict[dbTargetSeq]
                dDict[(dNewAlignData["GeneStart"],dNewAlignData["GeneStop"])]=dNewAlignData
                
                print(dOldAlignData)
                #print(dNewAlignData)
                print(dDict)
                
                #DEBUG
                #print("--------------Replace--------------")
                #print("Old alignment distance : {} \n{}\n{}\n{}".format(fOldDistance,
                    #dOldAlignData["GeneString"],dOldAlignData["AlignString"],dOldAlignData["ReadString"]))
                #print("New alignment distance : {} \n{}\n{}\n{}".format(fNewDistance,
                    #dNewAlignData["GeneString"],dNewAlignData["AlignString"],dNewAlignData["ReadString"]))
                #for sKey in dDict:
                    #if sKey!=(dNewAlignData["GeneStart"],dNewAlignData["GeneStop"]):
                        #print("Unchanged alignment distance : {} \n{}\n{}\n{}".format(dDict[sKey]["AlignIdentity"],
                        #dDict[sKey]["GeneString"],dDict[sKey]["AlignString"],dDict[sKey]["ReadString"]))
                #print("--------------/Replace--------------")
                #/DEBUG
                
                ##Remove the OldAlignData
                dToRemoveAlignData=self.get_tr_dict_data_alt()[sLineName][dbAlignGeneCoord]
                dbToRemoveReadCoord=(dToRemoveAlignData["ReadStart"],dToRemoveAlignData["ReadStop"])
                self.remove_key_from_tr_dict_data_alt(sLineName,dbAlignGeneCoord)
                self.remove_key_from_tr_dict_data(sLineName,dbToRemoveReadCoord)
                
                #1- Erase old Read data concerning CurrentBlock
                dbNeighboursMissing=self.erase_lineBlock(sCurrentBlock,iLineIndex)
                                
                ##Add the new data (unchanged and new)
                tKey=sorted(dDict.keys())
                for sKey in tKey:
                    dbGeneCoord=(dDict[sKey]["GeneStart"],dDict[sKey]["GeneStop"])
                    dbReadCoord=(dDict[sKey]["ReadStart"],dDict[sKey]["ReadStop"])
                    self.add_key_to_tr_dict_data_alt(sLineName,dbGeneCoord,dDict[sKey])
                    self.add_key_to_tr_dict_data(sLineName,dbReadCoord,dDict[sKey])
                    if sKey==tKey[0]:
                        self.update_alignmentMatrix(dDict[sKey],dbNeighboursMissing[0],None,iLineIndex)
                    elif sKey==tKey[-1]:
                        self.update_alignmentMatrix(dDict[sKey],None,dbNeighboursMissing[1],iLineIndex)
                    else:
                        self.update_alignmentMatrix(dDict[sKey],iUpstreamPenalty,iDownstreamPenalty,iLineIndex)
                    ##2- Update Read data concerning NewBlock
                    #self.update_alignmentMatrix(dNewAlignData,dbNeighboursMissing,iLineIndex)
                return True
                
        return False
    
    def add_key_to_tr_dict_data_alt(self,sLineName,dbGeneCoord,dDictContent):
        self.tr_dict_data_alt[sLineName][dbGeneCoord]=dDictContent
    
    def add_key_to_tr_dict_data(self,sLineName,dbReadCoord,dDictContent):
        self.tr_dict_data[sLineName][dbReadCoord]=dDictContent    
    
    def remove_key_from_tr_dict_data_alt(self,sLineName,dbGeneCoord):
        del self.tr_dict_data_alt[sLineName][dbGeneCoord]
    
    def remove_key_from_tr_dict_data(self,sLineName,dbReadCoord):
        del self.tr_dict_data[sLineName][dbReadCoord]
    
    def update_alignmentMatrix(self,dDict,oMissingUpstream,oMissingDownstream,iLineIndex):
        ##dDict as : {'ReadStop': 778, 'GeneSize': 20825, 'ReadStart': 747, 'ReadString': 'CAATCCGCCACTCGGATAAGTATGTCTGTCAT', 'GeneId': 'ENSMUSG00000000827', 'GeneStop': 14470, 'ReadSize': 1294, 'ReadId': 'ch173_read4429_template_pass_BYK_CB_ONT_1_FAF04998_A', 'GeneStart': 14440, 'AlignIdentity': 0.90625, 'AlignString': '| |||||||||||| ||||||||| |||||||', 'AlignSize': 32, 'AlignGap': 1, 'Strand': '-', 'GeneString': 'CCATCCGCCACTCG-ATAAGTATGCCTGTCAT'}
        self.update_lineBlock(iLineIndex,dDict["GeneStart"],dDict["GeneStop"],oMissingUpstream,oMissingDownstream)
        
    def update_lineBlock(self,iLineIndex,iStart,iStop,iPreviousStartValue=None,iNextStopValue=None):
        for iColIndex in range(iStart,iStop+1):
            self.get_matrix()[iLineIndex][iColIndex]=1
        if iPreviousStartValue is not None:
            if self.get_matrix()[iLineIndex][iStart-1]==0:
                self.get_matrix()[iLineIndex][iStart-1]=iPreviousStartValue
            else:
                exit("Error 2070 : iPreviousStartValue is not empty")
        if iNextStopValue is not None:
            if self.get_matrix()[iLineIndex][iStop+1]==0:
                self.get_matrix()[iLineIndex][iStop+1]=iNextStopValue
            else:
                exit("Error 2076 : iPreviousStartValue is not empty")
    
    def erase_lineBlock(self,sBlockName,iLineIndex):
        dbBlockCoord=self.get_blockCoord(sBlockName)
        iBlockStart=dbBlockCoord[0]
        iBlockStop=dbBlockCoord[1]
        iBlockPreviousStartCoord=iBlockStart-1
        iBlockNextStopCoord=iBlockStop+1
        
        for iColIndex in range(iBlockStart,iBlockStop+1):
            self.get_matrix()[iLineIndex][iColIndex]=0
        
        iPreviousValue=self.get_matrix()[iLineIndex][iBlockPreviousStartCoord]
        if iPreviousValue<0:
            self.get_matrix()[iLineIndex][iBlockPreviousStartCoord]=0
        else:
            iPreviousValue=None
        
        iNextValue=self.get_matrix()[iLineIndex][iBlockNextStopCoord]
        if iNextValue<0:
            self.get_matrix()[iLineIndex][iBlockNextStopCoord]=0
        else:
            iNextValue=None
            
        return (iPreviousValue,iNextValue)

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
                        elif iPreviousValue==-1:
                            ## Specific case : gap of size -1 can't have two equivalence
                            for iThisIndex in range(len(tData)):
                                #print(tBlockName[iThisIndex],tLineModel[iPreviousIndex-1])
                                if tBlockName[iThisIndex]==tLineModel[iPreviousIndex-1]:
                                    tData[iThisIndex+1]=str(iPreviousValue)
                                    break
                            iPreviousValue=tLineModel[iInLineIndex]
                            iPreviousIndex=iInLineIndex
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
                #print(iIndex,iAddIndex)
            #if self.get_globalVector(iIndex)=="i" and self.get_globalVector(iIndex-1)=="=":
                #iAddIndex=iIndex-1
                #print(iIndex,iAddIndex)
            if self.get_globalVector(iIndex)=="=" and self.get_globalVector(iIndex+1)!="=":
                iAddIndex=iIndex
                #print(iIndex,iAddIndex)
            #if self.get_globalVector(iIndex) not in ["i","="] or (self.get_globalVector(iIndex)=="i" and self.get_globalVector(iIndex-1)=="="):
            if iAddIndex is not None:
                if iLastIndex is None:
                    iLastIndex=iAddIndex
                    continue
                if iLastIndex+REGROUP_GLOBALVECTOR_VALUE>=iAddIndex:
                    #print(iAddIndex,self.get_globalVector(iAddIndex),"grouped")
                    if iCurrentGroup not in dGroup2Index:
                        dGroup2Index[iCurrentGroup]=[iLastIndex,iAddIndex]
                    else:
                        dGroup2Index[iCurrentGroup].append(iAddIndex)
                else:
                    #print(iAddIndex,self.get_globalVector(iAddIndex),"not grouped")
                    iCurrentGroup+=1
                iLastIndex=iAddIndex
        #print(dGroup2Index)
        
        ## Divide hybrid group
        bApplyCorrection=True
        while bApplyCorrection:
            bApplyCorrection=False
            for iGroupIdIndex in range(len(dGroup2Index)):
                iGroupId=sorted(dGroup2Index)[iGroupIdIndex]
                if len(set([self.get_globalVector(X) for X in dGroup2Index[iGroupId] if self.get_globalVector(X) in [">","<"]]))!=1:
                    #print("Hybrid Group : ",[self.get_globalVector(X) for X in dGroup2Index[iGroupId]],dGroup2Index[iGroupId])
                    bApplyCorrection=True
                    iCurrentGroup+=1
                    sTag=None
                    tNewGroup=[]
                    for iIndex in dGroup2Index[iGroupId][::-1]:
                        #print(iIndex,tNewGroup)
                        #print(self.get_globalVector(iIndex))
                        if self.get_globalVector(iIndex) not in ["<",">"]:
                            #tNewGroup.append(iIndex)
                            tNewGroup=[iIndex]+tNewGroup
                        elif sTag is None:
                            #tNewGroup.append(iIndex)
                            tNewGroup=[iIndex]+tNewGroup
                            sTag=self.get_globalVector(iIndex)
                        elif self.get_globalVector(iIndex)==sTag:
                            #tNewGroup.append(iIndex)
                            tNewGroup=[iIndex]+tNewGroup
                        else:
                            dGroup2Index[iCurrentGroup]=list(tNewGroup)
                            print("Divided into : ",[self.get_globalVector(X) for X in dGroup2Index[iCurrentGroup]],dGroup2Index[iCurrentGroup])
                            sTag=self.get_globalVector(iIndex)
                            tNewGroup=[iIndex]
                    iCurrentGroup+=1
                    dGroup2Index[iCurrentGroup]=list(tNewGroup)
                    #print("Divided into : ",[self.get_globalVector(X) for X in dGroup2Index[iCurrentGroup]],dGroup2Index[iCurrentGroup])
                    del dGroup2Index[iGroupId]
                    break
        
        ##DEBUG
        #for iGroupId in sorted(dGroup2Index):
            #print([self.get_globalVector(X) for X in dGroup2Index[iGroupId]],dGroup2Index[iGroupId])
        
        ##ASSIGN_GROUP_VALUE
        for iGroupId in sorted(dGroup2Index):
            #print("iGroupId",dGroup2Index[iGroupId])
            if len(set([self.get_globalVector(X) for X in dGroup2Index[iGroupId] if self.get_globalVector(X) in [">","<"]]))!=1:
                exit("ERROR 1726 : different type of variation\n{}\t{}".format([self.get_globalVector(X) for X in dGroup2Index[iGroupId]],dGroup2Index[iGroupId]))
            else:
                sTag=[self.get_globalVector(X) for X in dGroup2Index[iGroupId] if self.get_globalVector(X) in [">","<"]][0]
            if len([self.get_globalVector(X) for X in dGroup2Index[iGroupId] if self.get_globalVector(X) not in [">","<"]])==1:
                sBlockId=[self.get_globalVector(X) for X in dGroup2Index[iGroupId] if self.get_globalVector(X) not in [">","<"]][0]
            else:
                sBlockId=sTag
            #print("sBlockId",sBlockId,"sTag",sTag)
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
                            #tMatrix[iLineIndex][iColIndex]=0 # turn it to 0 now
                            #break
                            #print(">",iNegativeValue,iNegativeValueIndex)
                    #print("C>>>",tMatrix[iLineIndex][iStartColIndex-2:iStopColIndex+2])
                    for iColIndex in range(iStartColIndex-1,iStopColIndex+1): #same reason
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
                    for iColIndex in range(iStartColIndex,iStopColIndex+2): #+2, because negative signal is AFTER the grouping windows
                        if tMatrix[iLineIndex][iColIndex]<0:
                            bHaveNegativeValue=True
                            iNegativeValue=tMatrix[iLineIndex][iColIndex]
                            iNegativeValueIndex=iColIndex
                    #print("C'>>>",tMatrix[iLineIndex][iStartColIndex-2:iStopColIndex+2])
                    for iColIndex in range(iStartColIndex,iStopColIndex+2): #+2, same reason
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
            #print("~~~~~~~~~~~~~~~~~")
            #print(sTrId,dTr2Data[sTrId][sorted(dTr2Data[sTrId])[0]]["Strand"])
            tCurrentVector=list(tBaseVector)
            tDbKeyList=[]
            tConfidence=[]
            bOverlappingHit=False
            ## Step 3 : check colinearity
            for iCurrentDbIndex in range(len(dTr2Data[sTrId])-1):
                dbCurrentKey=sorted(dTr2Data[sTrId])[iCurrentDbIndex]
                sStrand=dTr2Data[sTrId][dbCurrentKey]["Strand"]
                iCurrentReadSize=dTr2Data[sTrId][dbCurrentKey]["ReadSize"]
                iCurrentGeneStart=dTr2Data[sTrId][dbCurrentKey]["GeneStart"]
                iCurrentGeneStop=dTr2Data[sTrId][dbCurrentKey]["GeneStop"]
                iCurrentReadStart=dTr2Data[sTrId][dbCurrentKey]["ReadStart"]
                iCurrentReadStop=dTr2Data[sTrId][dbCurrentKey]["ReadStop"]
                
                if iCurrentDbIndex==0:
                    tDbKeyList.append(dbCurrentKey)
                    tConfidence.append(0)
                
                for iAnotherDbIndex in range(iCurrentDbIndex+1,len(dTr2Data[sTrId])):
                    dbAnotherKey=sorted(dTr2Data[sTrId])[iAnotherDbIndex]
                    iAnotherGeneStart=dTr2Data[sTrId][dbAnotherKey]["GeneStart"]
                    iAnotherGeneStop=dTr2Data[sTrId][dbAnotherKey]["GeneStop"]
                    iAnotherReadStart=dTr2Data[sTrId][dbAnotherKey]["ReadStart"]
                    iAnotherReadStop=dTr2Data[sTrId][dbAnotherKey]["ReadStop"]
                    
                    if iCurrentDbIndex==0:
                        tDbKeyList.append(dbAnotherKey)
                        tConfidence.append(0)
                    
                    bStrongConfidence=False
                    if sStrand=="+":
                        if iCurrentReadStop<iAnotherReadStart:
                            ## iCurrentRead positionned before iAnotherRead
                            if iCurrentGeneStop<iAnotherGeneStart:
                                ## iCurrentGene positionned before iAnotherRead
                                bStrongConfidence=True
                        elif iCurrentReadStart>iAnotherReadStop:
                            ## iCurrentRead positionned after iAnotherRead
                            if iCurrentGeneStart>iAnotherGeneStop:
                                ## iCurrentGene positionned before iAnotherRead
                                bStrongConfidence=True
                        else:
                            bOverlappingHit=True
                            #exit("Error 2594 : iCurrentReadStart==iAnotherReadStop")
                    else:
                        if iCurrentReadStop<iAnotherReadStart:
                            ## iCurrentRead positionned before iAnotherRead
                            if iCurrentGeneStart>iAnotherGeneStop:
                                ## iCurrentGene positionned before iAnotherRead
                                bStrongConfidence=True
                        elif iCurrentReadStart>iAnotherReadStop:
                            ## iCurrentRead positionned after iAnotherRead
                            if iCurrentGeneStop<iAnotherGeneStart:
                                ## iCurrentGene positionned before iAnotherRead
                                bStrongConfidence=True
                        else:
                            bOverlappingHit=True
                            #exit("Error 2607 : iCurrentReadStart==iAnotherReadStop")
                        
                    if bStrongConfidence:
                        tConfidence[iCurrentDbIndex]+=1
                        tConfidence[iAnotherDbIndex]+=1
                        
            ## Step 2 : dismiss read where there is trouble into the confidence
            if max(tConfidence)!=min(tConfidence):
                if not bOverlappingHit:
                    print("Warning : read dismissed. Problem into the colinearity")
                else:
                    print("Warning : read dismissed. Overlaping hit")
                print("dbKey : {}".format(tDbKeyList))
                print("confidence : {}".format(tConfidence))
                continue
            
            ## Step 3 : fill the vector
            for iDbIndex in range(len(dTr2Data[sTrId])):
                dbCurrentKey=sorted(dTr2Data[sTrId])[iDbIndex]
                iCurrentReadSize=dTr2Data[sTrId][dbCurrentKey]["ReadSize"]
                iCurrentReadStart=dTr2Data[sTrId][dbCurrentKey]["ReadStart"]
                iCurrentReadStop=dTr2Data[sTrId][dbCurrentKey]["ReadStop"]
                iCurrentGeneSize=dTr2Data[sTrId][dbCurrentKey]["GeneSize"]
                iCurrentGeneStart=dTr2Data[sTrId][dbCurrentKey]["GeneStart"]
                iCurrentGeneStop=dTr2Data[sTrId][dbCurrentKey]["GeneStop"]
                sStrand=dTr2Data[sTrId][dbCurrentKey]["Strand"]
                
                #print("------------")
                #print("ReadSize",iCurrentReadSize,"GeneSize",iCurrentGeneSize)
                #print("dbCurrentKey",iCurrentReadStart,iCurrentReadStop,iCurrentGeneStart,iCurrentGeneStop)
                #print("TrData",dTr2Data[sTrId][dbCurrentKey])
                
                for iIndex in range(iCurrentGeneStart,iCurrentGeneStop+1):
                    if iIndex<iCurrentGeneSize:
                        tCurrentVector[iIndex]=1
                    else:
                        exit("Error 2257 : iIndex>=iCurrentGeneSize")
                if sum(tCurrentVector[iCurrentGeneStart:iCurrentGeneStop+1])!=len(tCurrentVector[iCurrentGeneStart:iCurrentGeneStop+1]):
                    exit("Error 3311 : Vector not correcting filled")
                        
                #print("Vector",tCurrentVector[iCurrentGeneStart-5:iCurrentGeneStop+1+5])
                #print("VectorCoord",iCurrentGeneStart-5,iCurrentGeneStop+1+5)
                
                bIsNotFirst=iDbIndex!=0
                if bIsNotFirst:
                    dbPreviousKey=sorted(dTr2Data[sTrId])[iDbIndex-1]
                    iPreviousReadSize=dTr2Data[sTrId][dbPreviousKey]["ReadSize"]
                    iPreviousGeneStart=dTr2Data[sTrId][dbPreviousKey]["GeneStart"]
                    iPreviousGeneStop=dTr2Data[sTrId][dbPreviousKey]["GeneStop"]
                    iPreviousReadStart=dTr2Data[sTrId][dbPreviousKey]["ReadStart"]
                    iPreviousReadStop=dTr2Data[sTrId][dbPreviousKey]["ReadStop"]
                else:
                    dbPreviousKey=None
                    iPreviousReadStart=0
                    iPreviousReadStop=0
                    iPreviousGeneStart=None
                    iPreviousGeneStop=None
                
                bIsNotLast=iDbIndex!=len(dTr2Data[sTrId])-1
                if bIsNotLast:
                    dbNextKey=sorted(dTr2Data[sTrId])[iDbIndex+1]
                    iNextReadSize=dTr2Data[sTrId][dbNextKey]["ReadSize"]
                    iNextGeneStart=dTr2Data[sTrId][dbNextKey]["GeneStart"]
                    iNextGeneStop=dTr2Data[sTrId][dbNextKey]["GeneStop"]
                    iNextReadStart=dTr2Data[sTrId][dbNextKey]["ReadStart"]
                    iNextReadStop=dTr2Data[sTrId][dbNextKey]["ReadStop"]
                else:
                    dbNextKey=None
                    iNextReadStart=iCurrentReadSize
                    iNextReadStop=iCurrentReadSize
                    iNextGeneStart=None
                    iNextGeneStop=None
                    
                bPreviousGap=False
                bNextGap=False
                
                #print("------------")
                #print("ReadSize",iCurrentReadSize,"sStrand",sStrand)
                #print("dbCurrentKey",iCurrentReadStart,iCurrentReadStop,iCurrentGeneStart,iCurrentGeneStop)
                #print("dbPreviousKey",iPreviousReadStart,iPreviousReadStop,iPreviousGeneStart,iPreviousGeneStop)
                #print("dbNextKey",iNextReadStart,iNextReadStop,iNextGeneStart,iNextGeneStop)
                                    
                if iPreviousReadStop+1==iCurrentReadStart:
                    ## Contigous alignment on the read, do nothing
                    pass
                elif iPreviousReadStop<iCurrentReadStart:
                    ## There is a gap
                    bPreviousGap=True
                    iPreviousGap=-(iCurrentReadStart-iPreviousReadStop-1)
                    #print("iPreviousGap",iPreviousGap)
                elif dbPreviousKey is None:
                    pass
                elif iPreviousReadStop>iCurrentReadStart:
                    ## Overlapping alignment on the Read. Not possible ?
                    exit("Error 2695 : overlapping alignment on the Read")
                else:
                    exit("FATAL 2697 : This line must never happend")
                
                if iCurrentReadStop+1==iNextReadStart:
                    ## Contigous alignment on the read, do nothing
                    pass
                elif iCurrentReadStop<iNextReadStart:
                    ## There is a gap
                    bNextGap=True
                    iNextGap=-(iNextReadStart-iCurrentReadStop-1)
                    #print("iNextGap",iNextGap)
                elif dbNextKey is None:
                    pass
                elif iCurrentReadStop>iNextReadStart:
                    ## Overlapping alignment on the Read. Not possible ?
                    exit("Error 2709 : overlapping alignment on the Read")
                else:
                    exit("FATAL 2710 : This line must never happend")
                
                if bPreviousGap: 
                    if sStrand=="+":
                        if iCurrentGeneStart-1>=0:
                            tCurrentVector[iCurrentGeneStart-1]=iPreviousGap
                        else:
                            tCurrentVector[0]=iPreviousGap+1
                        #iCurrentTargetIndex=max(iCurrentGeneStart-1,0)
                        ##print(iCurrentTargetIndex)
                        #tCurrentVector[iCurrentTargetIndex]=iPreviousGap
                        if dbPreviousKey is not None:
                            tCurrentVector[iPreviousGeneStop+1]=iPreviousGap
                    else:
                        if iCurrentGeneStop+1<=iCurrentGeneSize-1:
                            tCurrentVector[iCurrentGeneStop+1]=iPreviousGap
                        else:
                            tCurrentVector[iCurrentGeneSize-1]=iPreviousGap+1
                        #iCurrentTargetIndex=min(iCurrentGeneStop+1,iCurrentGeneSize-1)
                        ##print(iCurrentTargetIndex)
                        #tCurrentVector[iCurrentTargetIndex]=iPreviousGap
                        if dbPreviousKey is not None:
                            tCurrentVector[iPreviousGeneStart-1]=iPreviousGap
                
                if bNextGap: 
                    if sStrand=="+":
                        if iCurrentGeneStop+1<=iCurrentGeneSize-1:
                            tCurrentVector[iCurrentGeneStop+1]=iNextGap
                        else:
                            tCurrentVector[iCurrentGeneSize-1]=iNextGap+1
                        #iCurrentTargetIndex=min(iCurrentGeneStop+1,iCurrentGeneSize-1)
                        ##print(iCurrentTargetIndex)
                        #tCurrentVector[iCurrentTargetIndex]=iNextGap
                        if dbNextKey is not None:
                            tCurrentVector[iNextGeneStart-1]=iNextGap
                    else:
                        if iCurrentGeneStart-1>=0:
                            tCurrentVector[iCurrentGeneStart-1]=iNextGap
                        else:
                            tCurrentVector[0]=iNextGap+1
                        #iCurrentTargetIndex=max(iCurrentGeneStart-1,0)
                        ##print(iCurrentTargetIndex)
                        #tCurrentVector[iCurrentTargetIndex]=iNextGap
                        if dbNextKey is not None:
                            tCurrentVector[iNextGeneStop+1]=iNextGap
                
            tMatrix.append(tCurrentVector)
            self.set_matrix_lineName(sTrId,True)
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
    
    def parse_xmlFile(self,sPathFile,sGeneId):
        dRead2Coord={}
        dRead2GeneCoord={}
        for sLine in open(sPathFile):
            sStripedLine=sLine.strip()
            tLine=re.split("[></]",sStripedLine)
            tContentLine=[X for X in tLine if X!=""]
            if len(tContentLine)==1:
                continue
            sTagContent=tContentLine[0]
            sCoreContent=tContentLine[1]
            if sTagContent=="BlastOutput_query-def":
                sGeneId=sCoreContent
            if sTagContent=="BlastOutput_query-len":
                sGeneSize=sCoreContent
            if sTagContent=="Hit_id":
                sReadId=sCoreContent
                if sReadId not in dRead2Coord:
                    dRead2Coord[sReadId]={}
                    dRead2GeneCoord[sReadId]={}
            if sTagContent=="Hit_len":
                sReadSize=sCoreContent
            if sTagContent=="Hsp_query-from":
                sGene_HitStart=sCoreContent
            if sTagContent=="Hsp_query-to":
                sGene_HitStop=sCoreContent
            if sTagContent=="Hsp_hit-from":
                sRead_HitStart=sCoreContent
            if sTagContent=="Hsp_hit-to":
                sRead_HitStop=sCoreContent
            if sTagContent=="Hsp_hit-frame":
                #sRead_Strand=sCoreContent
                if sCoreContent=="1":
                    sRead_Strand="+"
                else:
                    sRead_Strand="-"
                    sRead_HitStart, sRead_HitStop = sRead_HitStop, sRead_HitStart
            if sTagContent=="Hsp_identity":
                sAlign_Identity=sCoreContent
            if sTagContent=="Hsp_gaps":
                sAlign_Gap=sCoreContent
            if sTagContent=="Hsp_align-len":
                sAlign_Length=sCoreContent
            if sTagContent=="Hsp_qseq":
                sGene_String=sCoreContent
            if sTagContent=="Hsp_hseq":
                sRead_String=sCoreContent
            if sTagContent=="Hsp_midline":
                sAlign_String=sCoreContent
                dRef2Content={
                    "ReadId":sReadId,
                    "ReadSize":int(sReadSize),
                    "ReadStart":int(sRead_HitStart),
                    "ReadStop":int(sRead_HitStop),
                    "Strand":sRead_Strand,
                    "GeneId":sGeneId,
                    "GeneSize":int(sGeneSize),
                    "GeneStart":int(sGene_HitStart),
                    "GeneStop":int(sGene_HitStop),
                    #"AlignIdentity":int(sAlign_Identity),
                    "AlignIdentity":len([X for X in range(len(sRead_String)) if sRead_String[X]==sGene_String[X]])/len(sRead_String),
                    "AlignGap":int(sAlign_Gap),
                    "AlignSize":int(sAlign_Length),
                    "GeneString":sGene_String,
                    "ReadString":sRead_String,
                    "AlignString":sAlign_String
                    }
                dRead2Coord[sReadId][(dRef2Content["ReadStart"],dRef2Content["ReadStop"])]=dRef2Content
                dRead2GeneCoord[sReadId][(min(dRef2Content["GeneStart"],dRef2Content["GeneStop"]),max(dRef2Content["GeneStart"],dRef2Content["GeneStop"]))]=dRef2Content
                
        self.set_size(int(dRef2Content["GeneSize"]))
        self.set_tr_dict_data(dRead2Coord)
        self.set_tr_dict_data_alt(dRead2GeneCoord)
    
    def parse_pafFile(self,sPathFile,sGeneId):
        dRead2Coord={}
        dRead2GeneCoord={}
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
                dRead2GeneCoord[sCurrentTrId]={}
            dRead2Coord[sCurrentTrId][(dRef2Content["ReadStart"],dRef2Content["ReadStop"])]=dRef2Content
            dRead2GeneCoord[sCurrentTrId][(min(dRef2Content["GeneStart"],dRef2Content["GeneStop"]),max(dRef2Content["GeneStart"],dRef2Content["GeneStop"]))]=dRef2Content
            
            #print(dRead2Coord)
            #print(dRead2GeneCoord)
            #exit()
            
        self.set_size(int(dRef2Content["GeneSize"]))
        self.set_tr_dict_data(dRead2Coord)
        self.set_tr_dict_data_alt(dRead2GeneCoord)

    def set_tr_dict_data_alt(self,dDict):
        self.tr_dict_data_alt=dDict
    
    def get_tr_dict_data_alt(self):
        return self.tr_dict_data_alt

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
