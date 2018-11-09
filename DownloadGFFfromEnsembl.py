# coding: utf-8
"""Python3.6"""

import os
import sys
from optparse import OptionParser
import urllib
import ast
import urllib.request
import json
import time

sCurrentVersionScript="v2"
########################################################################
'''
V2-2018/07/20
Explicit outputfile for fasta and gff

V1-2018/07/13
Download data for all transcript of a given GeneId and recreate the gff

python DownloadGFFfromEnsembl.py -t TARGETID -f FASTAFILE -g GFFFILE

TARGETID: gene of interest's EnsemblId
FASTAFILE: output path for fasta file
GFFFILE: output path for gff file
'''
########################################################################
#Options
parser = OptionParser(conflict_handler="resolve")
parser.add_option("-t","--targetid", dest="targetid")
parser.add_option("-f","--fastafile", dest="fastafile")
parser.add_option("-g","--gfffile", dest="gfffile")

(options, args) = parser.parse_args()

sTargetId=options.targetid
if not sTargetId:
    sys.exit("Error : no targetid -t defined, process broken")

sFastaFile=options.fastafile
if not sTargetId:
    sFastaFile=sTargetId+".fa"
    print("Warning : no fastafile -f defined, default : {}".format(sFastaFile))
    
sGffFile=options.gfffile
if not sTargetId:
    sGffFile=sTargetId+".gff"
    print("Warning : no gfffile -g defined, default : {}".format(sGffFile))


########################################################################
#CONSTANT
#DEBUG_GID="ENSMUSG00000027273"
VERBOSE=False

SERVER="http://rest.ensembl.org"
SERVER_ERROR429=429

GFF_FORMAT=("seqname","source","feature","start","end","score","strand","frame","attribute")

TRANSCRIPT_BIOTYPE=("protein_coding")
########################################################################
#CLASS
class EnsemblRestClient(object):
    """
    NOTE: methods __init__ and perform_rest_action has been taken from the Ensemble REST API website.
    """

    def __init__(self, sServer=SERVER):
        self.server = sServer
        
    def launch_request(self,sRequestString, bVerbose=VERBOSE):
        oData=None
        try:
            if bVerbose:
                print("Request:"+self.server+sRequestString)
            oRequest = urllib.request.Request(self.server + sRequestString)
            oResponse = urllib.request.urlopen(oRequest)
            oContent = oResponse.read()
            sContent=oContent.decode("utf-8")   
            try:         
                oData=json.loads(sContent)
            except json.decoder.JSONDecodeError:
                oData=sContent

        except urllib.request.HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == SERVER_ERROR429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    oData=self.launch_request(sRequestString)
            else:
                sys.stderr.write(
                    'Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))    
        return oData
            
    def get_geneData(self,sGeneId):
        sRequestSuffix="/lookup/id/"+sGeneId
        sRequestSuffix+="?content-type=application/json" #results as json format
        sRequestSuffix+=";expand=1" #expand all item linked to sGeneId (transcript, exon, etc.)
        sRequestSuffix+=";utr=1" #Get the utr data
        dGene=self.launch_request(sRequestSuffix)
        return dGene
        
    def get_geneSequence(self,sGeneId):
        sRequestSuffix="/sequence/id/"+sGeneId
        sRequestSuffix+="?content-type=text/plain" #results as json format
        sGeneSeq=self.launch_request(sRequestSuffix)
        return sGeneSeq

class EnsemblGene():
    
    def __init__(self,dDict):
        dbData=get_itemInfo(dDict)
        self.set_seqname(dbData[0])
        self.set_source(dbData[1])
        self.set_feature(dbData[2])
        self.set_start(dbData[3])
        self.set_end(dbData[4])
        self.set_score(dbData[5])
        self.set_strand(dbData[6])
        self.set_frame(dbData[7])
        self.set_attribute(dbData[8])
        self.set_emptyTranscriptList()
        
        for dTr in dDict["Transcript"]:
            if dTr["biotype"] in TRANSCRIPT_BIOTYPE:
                oTranscript=EnsemblTranscript(dTr,self)
                self.add_transcript(oTranscript)
                
    def describe(self):
        sContent=""
        for sGffPart in GFF_FORMAT:
            if sContent!="":
                sContent+="\t"
            sContent+=str(getattr(self,sGffPart))
        sContent+="\n"
        for oTranscript in sorted(self.get_transcriptList(),key=lambda x: x.get_start()):
            sContent+=oTranscript.describe()

        return sContent
    
    def set_emptyTranscriptList(self):
        self.transcriptList=[]
    
    def add_transcript(self,oObject):
        self.transcriptList.append(oObject)
        
    def get_transcriptList(self):
        return self.transcriptList
    
    def set_seqname(self,sString):
        self.seqname=sString
    
    def set_source(self,sString):
        self.source=sString
        
    def set_feature(self,sString):
        self.feature=sString
    
    def set_start(self,iValue):
        self.start=iValue
    
    def set_end(self,iValue):
        self.end=iValue
        
    def set_score(self,sString):
        self.score=sString
        
    def set_strand(self,sString):
        self.strand=sString
        
    def set_strand(self,sString):
        self.strand=sString
        
    def set_frame(self,sString):
        self.frame=sString
        
    def set_attribute(self,sString):
        self.attribute=sString
        
    def get_seqname(self):
        return self.seqname
    
    def get_source(self):
        return self.source
        
    def get_feature(self):
        return self.feature
    
    def get_start(self):
        return self.start
    
    def get_end(self):
        return self.end
        
    def get_score(self):
        return self.score
        
    def get_strand(self):
        return self.strand
        
    def get_strand(self):
        return self.strand
        
    def get_frame(self):
        return self.frame
        
    def get_attribute(self):
        return self.attribute
        
class EnsemblTranscript:
    
    def __init__(self,dDict,oParent):
        dbData=get_itemInfo(dDict)
        self.set_seqname(dbData[0])
        self.set_source(dbData[1])
        self.set_feature(dbData[2])
        self.set_start(dbData[3])
        self.set_end(dbData[4])
        self.set_score(dbData[5])
        self.set_strand(dbData[6])
        self.set_frame(dbData[7])
        self.set_attribute(dbData[8])
        self.add_attribute(oParent.get_attribute())
        self.set_emptyExonList()
        
        for dExon in dDict["Exon"]:
            oExon=EnsemblExon(dExon,self)
            self.add_exon(oExon)
    
                    
    def describe(self):
        sContent=""
        for sGffPart in GFF_FORMAT:
            if sContent!="":
                sContent+="\t"
            sContent+=str(getattr(self,sGffPart))
        sContent+="\n"
        for oExon in sorted(self.get_exonList(),key=lambda x: x.get_start()):
            sContent+=oExon.describe()

        return sContent
    
    def set_emptyExonList(self):
        self.exonList=[]
    
    def add_exon(self,oObject):
        self.exonList.append(oObject)
        
    def get_exonList(self):
        return self.exonList
    
    def set_seqname(self,sString):
        self.seqname=sString
    
    def set_source(self,sString):
        self.source=sString
        
    def set_feature(self,sString):
        self.feature=sString
    
    def set_start(self,iValue):
        self.start=iValue
    
    def set_end(self,iValue):
        self.end=iValue
        
    def set_score(self,sString):
        self.score=sString
        
    def set_strand(self,sString):
        self.strand=sString
        
    def set_strand(self,sString):
        self.strand=sString
        
    def set_frame(self,sString):
        self.frame=sString
        
    def set_attribute(self,sString):
        self.attribute=sString
        
    def add_attribute(self,sString):
        self.attribute+=";"+sString
        
    def get_seqname(self):
        return self.seqname
    
    def get_source(self):
        return self.source
        
    def get_feature(self):
        return self.feature
    
    def get_start(self):
        return self.start
    
    def get_end(self):
        return self.end
        
    def get_score(self):
        return self.score
        
    def get_strand(self):
        return self.strand
        
    def get_strand(self):
        return self.strand
        
    def get_frame(self):
        return self.frame
        
    def get_attribute(self):
        return self.attribute
            
class EnsemblExon:
    
    def __init__(self,dDict,oParent):
        dbData=get_itemInfo(dDict)
        self.set_seqname(dbData[0])
        self.set_source(oParent.get_source())
        self.set_feature(dbData[2])
        self.set_start(dbData[3])
        self.set_end(dbData[4])
        self.set_score(dbData[5])
        self.set_strand(dbData[6])
        self.set_frame(dbData[7])
        self.set_attribute(dbData[8])
        self.add_attribute(oParent.get_attribute())
    
    def describe(self):
        sContent=""
        for sGffPart in GFF_FORMAT:
            if sContent!="":
                sContent+="\t"
            sContent+=str(getattr(self,sGffPart))
        sContent+="\n"

        return sContent
    
    def set_seqname(self,sString):
        self.seqname=sString
    
    def set_source(self,sString):
        self.source=sString
        
    def set_feature(self,sString):
        self.feature=sString
    
    def set_start(self,iValue):
        self.start=iValue
    
    def set_end(self,iValue):
        self.end=iValue
        
    def set_score(self,sString):
        self.score=sString
        
    def set_strand(self,sString):
        self.strand=sString
        
    def set_strand(self,sString):
        self.strand=sString
        
    def set_frame(self,sString):
        self.frame=sString
        
    def set_attribute(self,sString):
        self.attribute=sString
        
    def add_attribute(self,sString):
        self.attribute+=";"+sString
        
    def get_seqname(self):
        return self.seqname
    
    def get_source(self):
        return self.source
        
    def get_feature(self):
        return self.feature
    
    def get_start(self):
        return self.start
    
    def get_end(self):
        return self.end
        
    def get_score(self):
        return self.score
        
    def get_strand(self):
        return self.strand
        
    def get_strand(self):
        return self.strand
        
    def get_frame(self):
        return self.frame
        
    def get_attribute(self):
        return self.attribute

        
########################################################################
#Function 
def get_itemInfo(dDict,oParent=None):
    sStrand="+"
    if dDict["strand"]==-1:
        sStrand="-"
    sObjectType=dDict["object_type"]
    sObjectTypeId=sObjectType.lower()
    sObjectTypeId=sObjectTypeId.replace("_prime","'")
    sObjectTypeId=sObjectTypeId.replace("three","3")
    sObjectTypeId=sObjectTypeId.replace("five","5")
    if "UTR" in sObjectType:
        sObjectType="UTR"
    try:
        sSource=dDict["source"]
    except KeyError:
        sSource="."
    return(dDict["seq_region_name"],
            sSource,
            sObjectType,
            dDict["start"],
            dDict["end"],
            ".",
            sStrand,
            ".",
            "{}_id {}".format(sObjectTypeId,dDict["id"])
            )

def WriteGFF(sFile,oObject):
    sContent=oObject.describe()
    WriteFile(sFile,sContent)
    
def WriteMonoFasta(sSeqTitle,sFileName,sSeqContent):
    sContent=">"+sSeqTitle+"\n"
    iBase=1
    iLineSize=60
    while iBase*iLineSize<len(sSeqContent):
        sContent+=sSeqContent[(iBase-1)*iLineSize:iBase*iLineSize]+"\n"
        iBase+=1
    sContent+=sSeqContent[(iBase-1)*iLineSize:]+"\n"
    WriteFile(sFileName,sContent)

def WriteFile(sFile,sContent,bVerbose=VERBOSE):
    FILE=open(sFile,"w")
    FILE.write(sContent)
    FILE.close()
    if bVerbose:
        print("File writed : "+sFile)
    
########################################################################
#MAIN
if __name__ == "__main__":
    er = EnsemblRestClient()
    dGeneData=er.get_geneData(sTargetId)
    sGeneSeq=er.get_geneSequence(sTargetId)
    oEnsemblGene=EnsemblGene(dGeneData)
    WriteGFF(sGffFile,oEnsemblGene)
    WriteMonoFasta(sTargetId,sFastaFile,sGeneSeq)
    #print("DONE")
