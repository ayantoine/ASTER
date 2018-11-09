# coding: utf-8
"""Python3.6"""

import os
import sys
from optparse import OptionParser
#import matplotlib.pyplot as plt

sCurrentVersionScript="v5"
########################################################################
'''
V5-2018/10/24
Added:Threshold for ReadCover and GapedGeneCover (expressed in %)

V4-2018/10/18
Adapt to new ResultFile
V3-2018/10/10
If no targetRead, take only TrId in sourceFile if column Index3 is not -1 (%ReadCover, set at -1 if chimeric)
V2-2018/09/25
Extract targetread or read them from ExtractMergeLine.py output
'''
########################################################################
#Options
parser = OptionParser(conflict_handler="resolve")
parser.add_option("-f","--fastafile", dest="fastafile")
parser.add_option("-t","--targetread", dest="targetread")
parser.add_option("-s","--sourcefile", dest="sourcefile")
parser.add_option("-o","--outputfile", dest="outputfile")
parser.add_option("-r","--readcoverthreshold", dest="readcoverthreshold")
parser.add_option("-g","--gapedgenecoverthreshold", dest="gapedgenecoverthreshold")

(options, args) = parser.parse_args()

sFastaFile=options.fastafile
if not sFastaFile:
    sys.exit("Error : no fastafile -f defined, process broken")
    
sTargetRead=options.targetread
sSourceFile=options.sourcefile
if not sTargetRead and not sSourceFile:
    sys.exit("Error : no targetread -t or sourcefile -s defined, process broken")
elif sTargetRead and sSourceFile:
    sys.exit("Error : only one among targetread -t and sourcefile -s must be defined, process broken")
elif sTargetRead:
    bTargetRead=True
elif sSourceFile:
    bTargetRead=False
else:
    sys.exit("FATAL : L41")

sOutputFile=options.outputfile
if not sOutputFile:
    sys.exit("Error : no outputfile -o defined, process broken")

sReadCoverThreshold=options.readcoverthreshold
try:
    fReadCoverThreshold=float(sReadCoverThreshold)
except:
    fReadCoverThreshold=0.0
    print("Warning :  no readcoverthreshold -r defined, default value will be {} %")

sGeneCoverThreshold=options.gapedgenecoverthreshold
try:
    fGeneCoverThreshold=float(sGeneCoverThreshold)
except:
    fGeneCoverThreshold=0.0
    print("Warning :  no gapedgenecoverthreshold -g defined, default value will be {} %")


########################################################################
#DEBUG
#print("sFastaFile",sFastaFile)
#print("sTargetRead",sTargetRead)
#print("sSourceFile",sSourceFile)
#print("sOutputFile",sOutputFile)

########################################################################
#CONSTANT
COLUMN_READCOVER_INDEX=5
COLUMN_GAPEDGENECOVER_INDEX=8
COLUMN_READID_INDEX=3
#BAD_READCOVER_VALUE=-1

READCOVER_THRESHOLD=fReadCoverThreshold
GENECOVER_THRESHOLD=fGeneCoverThreshold
########################################################################
#Function
def LoadFastaFile(sPath):
    dDict={}
    sSeqName=""
    sSeqContent=""
    for sLine in open(sPath):
        sContent=sLine.strip()
        if sContent[0]==">":
            if sSeqName!="":
                try:
                    dDict[sSeqName].append(sSeqContent)
                except KeyError:
                    dDict[sSeqName]=[sSeqContent]
            sSeqName=sContent.replace(">","")
            sSeqContent=""
        else:
            sSeqContent+=sContent
    try:
        dDict[sSeqName].append(sSeqContent)
    except KeyError:
        dDict[sSeqName]=[sSeqContent]
    return dDict

def GetTrIdListFromSourceFile(sFile):
    tList=[]
    bHeader=True
    for sLine in open(sFile):
        if bHeader:
            bHeader=False
            continue
        tLine=sLine.split("\t")
        #print(tLine)
        #if float(tLine[COLUMN_READCOVER_INDEX])!=BAD_READCOVER_VALUE:
        if float(tLine[COLUMN_READCOVER_INDEX])>=READCOVER_THRESHOLD \
        and float(tLine[COLUMN_GAPEDGENECOVER_INDEX])>=GENECOVER_THRESHOLD:
            tList.append(tLine[COLUMN_READID_INDEX])
    return tList

def WriteFasta(dDict,sPath,tList=None):
    if tList!=None:
        tTarget=list(tList)
    else:
        tTarget=dDict.keys()
    FILE=open(sPath,"w")
    for sKey in tTarget:
        try:
            oContent=dDict[sKey]
        except KeyError:
            print("WARNING : target {} not present into the fasta".format(sKey))
        if len(oContent)!=1:
            print("WARNING : target {} is referenced {} time into the fasta".format(sKey,len(oContent)))
            continue
        FILE.write(">{}\n".format(sKey))
        FILE.write("{}\n".format(dDict[sKey][0]))
    FILE.close()

########################################################################
#MAIN
dFasta=LoadFastaFile(sFastaFile)
#print(len(dFasta))
#print(dFasta)
if bTargetRead:
    tTargetRead=sTargetRead.split(";")
else:
    tTargetRead=GetTrIdListFromSourceFile(sSourceFile)
WriteFasta(dFasta,sOutputFile,tTargetRead)

