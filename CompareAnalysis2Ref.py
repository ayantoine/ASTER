# coding: utf-8
"""Python3.6"""

import os
import sys
from optparse import OptionParser

from CurrentWorkflowClass import *

sCurrentVersionScript="v8"
########################################################################
'''
V8-2018/12/20
Rework all the logic. Centralize on megablast and exonerate

python FilterBlastResult.py -b BLASTFILE -t TARGETID -g GFFFILE 
-e EXONFILE -o OUTPUTFILE

BLASTFILE: path to the blast xml file result
TARGETID: gene of interest's EnsemblId
GFFFILE: path to the gff file of reference
EXONFILE: path to store the list of exon and all transcript associated
OUTPUTFILE: path to store the result matrix
'''
########################################################################
#Options
parser = OptionParser(conflict_handler="resolve")
parser.add_option("-b","--blastfile", dest="blastfile")
parser.add_option("-t","--targetid", dest="targetid")
parser.add_option("-g","--gfffile", dest="gfffile")
parser.add_option("-o","--outputfile", dest="outputfile")
parser.add_option("-e","--exonlistfile", dest="exonlistfile")

(options, args) = parser.parse_args()

sBlastFile=options.blastfile
if not sBlastFile:
    sys.exit("Error : no blastfile -b defined, process broken")

sOutputFile=options.outputfile
if not sOutputFile:
    sOutputFile="Result_"+os.path.basename(sInputFile).split(".")[0]+".tsv"
    print("Warning : no outputfile -o defined, default {}".format(sOutputFile))

sExonListFile=options.exonlistfile
if not sExonListFile:
    sExonListFile="ExonList_"+os.path.basename(sInputFile).split(".")[0]+".txt"
    print("Warning : no exonlistfile -e defined, default {}".format(sExonListFile))

sTargetId=options.targetid
if not sTargetId:
    sys.exit("Error : no targetid -t defined, process broken")
    
sGffFile=options.gfffile
if not sGffFile:
    sys.exit("Error : no gfffile -g defined, process broken")

########################################################################
#Function 
def WriteFile(sPath,sContent):
    FILE=open(sPath,"w")
    FILE.write(sContent)
    FILE.close()
  
#def ExtractRead2Strand(sPath):
    #dResult={}
    #bHeader=True
    #for sLine in open(sPath):
        #if bHeader:
            #bHeader=False
            #continue
        #tLine=sLine.split("\t")
        #iGeneStrand=int(tLine[1])
        #iReadStrand=int(tLine[2])
        #sReadId=tLine[3]
        #dResult[sReadId]={"geneStrand":iGeneStrand,"readStrand":iReadStrand}
    #return dResult
    
########################################################################
#MAIN
if __name__ == "__main__":
    oRefContent=TranscriptContent(sTargetId,sGffFile)
    print("Ref contains {} exons along {} isoformes".format(len(oRefContent.get_exonList()),len(oRefContent.get_transcript())))
    sExonList=oRefContent.describe_exon()
    WriteFile(sExonListFile,sExonList)
    
    sTool="Blast"
    oFileContent=BlastContent(sBlastFile)
    #oFileContent.purge_alignList(["ch27_read13167_template_pass_BYK_CB_ONT_1_FAF04998_A"])
    #oFileContent.describe_AllHsp()
    print("{} length : {}nt".format(oFileContent.query_id,oFileContent.query_size))
    print("{} associated reads : {}".format(sTool,len(oFileContent.align_list)))
        
    ##TSV exon covering file
    oFileContent.compute_self_cover_alignments()
    oFileContent.define_position_alignments()
    oFileContent.assign_exonCovering(oRefContent)
    sContent=oFileContent.describe_alignment_header(oRefContent,False)
    sContent+=oFileContent.describe_alignment(oRefContent,False)
    WriteFile(sOutputFile,sContent)
    
    ###TSV Ref2Read coord
    WriteFile(".".join(sOutputFile.split(".")[:-1])+".paf",oFileContent.export_data())
