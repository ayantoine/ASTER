# coding: utf-8
"""Python3.6"""

import os
import sys
import re
from optparse import OptionParser

#from CurrentWorkflowClass import *

sCurrentVersionScript="v1"
########################################################################
'''
V1-2018/12/20
Parse exonerate output. provide blast-xml-like output.

python ParseExonerateOutput.py -e EXONERATEFILE -b BLASTFILE -o OUTPUTFILE
-g GENEFILE -r READFILE

EXONERATEFILE: path to the exonerate txt file result
BLASTFILE: path to the blast tsv table
TARGETID: gene of interest's EnsemblId
OUTPUTFILE: path to store the result matrix
GENEFILE: path to the gene fasta file
READFILE: path to the read fasta file
'''
########################################################################
#Options
parser = OptionParser(conflict_handler="resolve")
parser.add_option("-e","--exoneratefile", dest="exoneratefile")
parser.add_option("-b","--blastfile", dest="blastfile")
parser.add_option("-g","--genefile", dest="genefile")
parser.add_option("-r","--readfile", dest="readfile")
parser.add_option("-o","--outputfile", dest="outputfile")

(options, args) = parser.parse_args()

sBlastFile=options.blastfile
if not sBlastFile:
    sys.exit("Error : no blastfile -b defined, process broken")

sExonerateFile=options.exoneratefile
if not sExonerateFile:
    sys.exit("Error : no exoneratefile -e defined, process broken")

sOutputFile=options.outputfile
if not sOutputFile:
    sOutputFile="Result_"+os.path.basename(sInputFile).split(".")[0]+".tsv"
    print("Warning : no outputfile -o defined, default {}".format(sOutputFile))

sGeneFile=options.genefile
if not sGeneFile:
    sys.exit("Error : no genefile -g defined, process broken")

sReadFile=options.readfile
if not sReadFile:
    sys.exit("Error : no readfile -r defined, process broken")

########################################################################
#CONSTANT
REVERSE_COMP={"A":"T","C":"G","G":"C","T":"A",
                "a":"t","c":"g","g":"c","t":"a",
                "N":"N","n":"n",
                ">":"<","<":">",".":".","-":"-"," ":" "
                }
RAWSCORE_THRESHOLD=150
########################################################################
#Function 
def WriteFile(sPath,sContent):
    FILE=open(sPath,"w")
    FILE.write(sContent)
    FILE.close()

def ExtractRead2Strand(sPath):
    dResult={}
    bHeader=True
    for sLine in open(sPath):
        if bHeader:
            bHeader=False
            continue
        tLine=sLine.split("\t")
        iGeneStrand=int(tLine[1])
        iReadStrand=int(tLine[2])
        sReadId=tLine[3]
        dResult[sReadId]={"geneStrand":iGeneStrand,"readStrand":iReadStrand}
    return dResult

def exonerate_reverseComp(sString):
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
    
    
def exonerate_reverseAlign(sString):
    tString=re.split("\+\+|-\+|\+-|--",sString)
    tResult=[]
    for iIndex in range(len(tString)-1,-1,-1):
        if "|" in tString[iIndex]:
            tResult.append(tString[iIndex][::-1])
        else:
            tResult.append(tString[iIndex])
    sResult="++".join(tResult)
    return sResult

def parse_exonerateFile(sFile,dGene2Fasta,sReadFastaFile,dRead2Strand):
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
                    
                    #print("dRead2Strand",dRead2Strand[sReadId]["geneStrand"],dRead2Strand[sReadId]["readStrand"])
                    #print("Exonerate",iGeneStrand,iReadStrand)
                    
                    if iGeneStrand==1:
                        if dRead2Strand[sReadId]["readStrand"]==iReadStrand:
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
                            print("Reject alignment from {} : strandRead is diverging from megablast analyze")
                    else:
                        if dRead2Strand[sReadId]["readStrand"]==-iReadStrand:
                            dResult[(sGeneId,iGeneSize)][(sReadId,iReadSize)][len(dResult[(sGeneId,iGeneSize)][(sReadId,iReadSize)])]={
                                        "frame":(-iGeneStrand,-iReadStrand),
                                        "query":exonerate_reverseComp(sGeneAlign),
                                        "query_start":iGeneStop-1,
                                        "query_end":iGeneStart-1,
                                        "sbjct":exonerate_reverseComp(sReadAlign),
                                        "sbjct_start":iReadStop-1,
                                        "sbjct_end":iReadStart-1,
                                        "align_data":exonerate_reverseAlign(sDataAlign),
                                        "raw_score":iRawScore
                                        }
                        else:
                            print("Reject alignment from {} : strandRead is diverging from megablast analyze")
                            
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

def refine_exonerateParsing(dDict):
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
                                
    return dDict

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
    
def xml_header(dDict):
    for dbKey in dDict:
        return """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.5.0+</BlastOutput_version>
  <BlastOutput_reference>Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), &quot;A greedy algorithm for aligning DNA sequences&quot;, J Comput Biol 2000; 7(1-2):203-14.</BlastOutput_reference>
  <BlastOutput_db></BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>"""+dbKey[0]+"""</BlastOutput_query-def>
  <BlastOutput_query-len>"""+str(dbKey[1])+"""</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>0</Parameters_expect>
      <Parameters_sc-match>0</Parameters_sc-match>
      <Parameters_sc-mismatch>0</Parameters_sc-mismatch>
      <Parameters_gap-open>0</Parameters_gap-open>
      <Parameters_gap-extend>0</Parameters_gap-extend>
      <Parameters_filter>-</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
<BlastOutput_iterations>
<Iteration>
  <Iteration_iter-num>1</Iteration_iter-num>
  <Iteration_query-ID>Query_1</Iteration_query-ID>
  <Iteration_query-def>"""+dbKey[0]+"""</Iteration_query-def>
  <Iteration_query-len>"""+str(dbKey[1])+"""</Iteration_query-len>
<Iteration_hits>
<Hit>
"""

def xml_core(dDict,FILE):
    for dbGeneKey in dDict:
        iKeyCounter=0
        for dbKey in dDict[dbGeneKey]:
            iKeyCounter+=1
            FILE.write("""  <Hit_num>"""+str(iKeyCounter)+"""</Hit_num>
  <Hit_id>"""+dbKey[0]+"""</Hit_id>
  <Hit_def></Hit_def>
  <Hit_accession>"""+dbKey[0]+"""</Hit_accession>
  <Hit_len>"""+str(dbKey[1])+"""</Hit_len>
  <Hit_hsps>
""")
            for iHspIndex in dDict[dbGeneKey][dbKey]:
                FILE.write("""    <Hsp>
      <Hsp_num>"""+str(iHspIndex)+"""</Hsp_num>
      <Hsp_bit-score>0</Hsp_bit-score>
      <Hsp_score>0</Hsp_score>
      <Hsp_evalue>0</Hsp_evalue>
      <Hsp_query-from>"""+str(dDict[dbGeneKey][dbKey][iHspIndex]["query_start"])+"""</Hsp_query-from>
      <Hsp_query-to>"""+str(dDict[dbGeneKey][dbKey][iHspIndex]["query_end"])+"""</Hsp_query-to>
      <Hsp_hit-from>"""+str(dDict[dbGeneKey][dbKey][iHspIndex]["sbjct_start"])+"""</Hsp_hit-from>
      <Hsp_hit-to>"""+str(dDict[dbGeneKey][dbKey][iHspIndex]["sbjct_end"])+"""</Hsp_hit-to>
      <Hsp_query-frame>"""+str(dDict[dbGeneKey][dbKey][iHspIndex]["frame"][0])+"""</Hsp_query-frame>
      <Hsp_hit-frame>"""+str(dDict[dbGeneKey][dbKey][iHspIndex]["frame"][-1])+"""</Hsp_hit-frame>
      <Hsp_identity>"""+str(dDict[dbGeneKey][dbKey][iHspIndex]["align_data"].count("|"))+"""</Hsp_identity>
      <Hsp_positive>"""+str(dDict[dbGeneKey][dbKey][iHspIndex]["align_data"].count("|"))+"""</Hsp_positive>
      <Hsp_gaps>"""+str(dDict[dbGeneKey][dbKey][iHspIndex]["sbjct"].count("-")+dDict[dbGeneKey][dbKey][iHspIndex]["query"].count("-"))+"""</Hsp_gaps>
      <Hsp_align-len>"""+str(len(dDict[dbGeneKey][dbKey][iHspIndex]["query"]))+"""</Hsp_align-len>
      <Hsp_qseq>"""+dDict[dbGeneKey][dbKey][iHspIndex]["query"]+"""</Hsp_qseq>
      <Hsp_hseq>"""+dDict[dbGeneKey][dbKey][iHspIndex]["sbjct"]+"""</Hsp_hseq>
      <Hsp_midline>"""+dDict[dbGeneKey][dbKey][iHspIndex]["align_data"]+"""</Hsp_midline>
    </Hsp>
""")
            FILE.write("""  </Hit_hsps>\n""")
    
def xml_tail():
    return """</Hit>
</Iteration_hits>
  <Iteration_stat>
    <Statistics>
      <Statistics_db-num>0</Statistics_db-num>
      <Statistics_db-len>0</Statistics_db-len>
      <Statistics_hsp-len>0</Statistics_hsp-len>
      <Statistics_eff-space>0</Statistics_eff-space>
      <Statistics_kappa>0</Statistics_kappa>
      <Statistics_lambda>0</Statistics_lambda>
      <Statistics_entropy>0</Statistics_entropy>
    </Statistics>
  </Iteration_stat>
</Iteration>
</BlastOutput_iterations>
</BlastOutput>"""
    
########################################################################
#MAIN
if __name__ == "__main__":
    ##Get major strand determined by megablast
    dRead2Strand=ExtractRead2Strand(sBlastFile)
    
    dGene2Fasta=Fasta2Dict(sGeneFile)
    if len(dGene2Fasta)!=1:
        exit("ERROR 383 : Many Ref in Ref fasta file")
    
    dRead2ExonerateData=parse_exonerateFile(sExonerateFile,dGene2Fasta,sReadFile,dRead2Strand)
    dRead2ExonerateData=refine_exonerateParsing(dRead2ExonerateData)
    #print(len(dRead2ExonerateData))
    #print(dRead2ExonerateData.keys())
    
    FILE=open(sOutputFile,"w")
    FILE.write(xml_header(dRead2ExonerateData))
    xml_core(dRead2ExonerateData,FILE)
    FILE.write(xml_tail())
    FILE.close()
