# coding: utf-8
"""Python3.6"""

import os
import sys
from optparse import OptionParser

from CurrentWorkflowClass import *

sCurrentVersionScript="v7"
########################################################################
'''
V7-2018/10/18
Added: Exonerate file can use X.megablast.tsv results to reject Hsp

V6-2018/09/28
Added: Exonerate file
V5-2018/09/20
Added: Externalise Class Object into separate file. Easiest to share with other script
V4-2018/09/03
Some hit on exon was considered as non-exon hit
V3-2018/08/23
oHsp end use specific attribut. More confidence than calcul
Merge with CompareYassRaw2Ref
#patch-2018/08/30
TranscriptContent store Gene coord/strand
correct_exonCoord use now the real Gene start
V2-2018/07/11
Use gff file for the Ref
V1-2018/06/04
Parse a xml blast result file and compare to exon list in target

python FilterBlastResult.py -b BLASTFILE -t TARGETID -g GFFFILE 
(-x EXONERATEFILE -y YASSFILE -e EXONFILE -o OUTPUTFILE -s MEGABLASTSTRAND) 

BLASTFILE: path to the blast xml file result
TARGETID: gene of interest's EnsemblId
GFFFILE: path to the gff file of reference
YASSFILE: path to the yass raw file result (alternativ to BLASTFILE)
EXONERATE: path to the exonerate raw file result (alternativ to BLASTFILE)
EXONFILE: path to store the list of exon and all transcript associated
OUTPUTFILE: path to store the result matrix
MEGABLASTSTRAND: path to the refined megablast result file to specify Hsp strand filtering
'''
########################################################################
#Options
parser = OptionParser(conflict_handler="resolve")
parser.add_option("-b","--blastfile", dest="blastfile")
parser.add_option("-y","--yassfile", dest="yassfile")
parser.add_option("-x","--exoneratefile", dest="exoneratefile")
parser.add_option("-r","--reffile", dest="reffile")
parser.add_option("-f","--fastafile", dest="fastafile")
parser.add_option("-t","--targetid", dest="targetid")
parser.add_option("-g","--gfffile", dest="gfffile")
parser.add_option("-o","--outputfile", dest="outputfile")
parser.add_option("-e","--exonlistfile", dest="exonlistfile")
parser.add_option("-s","--megablaststrand", dest="megablaststrand")

(options, args) = parser.parse_args()

sBlastFile=options.blastfile
sYassFile=options.yassfile
sExonerateFile=options.exoneratefile
sRefFile=options.reffile
sFastaFile=options.fastafile
tInputFile=[sBlastFile is not None, sYassFile is not None, sExonerateFile is not None]
if sum(tInputFile)==0:
    sys.exit("Error : no inputfile -b/-y/-x defined, process broken")
elif sum(tInputFile)!=1:
    sys.exit("Error : multiple inputfile -b/-y/-x defined, only one is allowed, process broken")
elif sBlastFile:
    sInputFile=sBlastFile
    #bBlastInput=True
elif sYassFile:
    sInputFile=sYassFile
    #bBlastInput=False
elif sExonerateFile:
    sInputFile=sExonerateFile
    if not sFastaFile:
        exit("Error : no fastafile -f defined. ReadFasta needed for exonerate parsing")
    if not sRefFile:
        exit("Error : no reffile -r defined. GeneFasta needed for exonerate parsing")
iTool=tInputFile.index(True) #0:Blast, 1:Yass, 2:Exonerate

sOutputFile=options.outputfile
if not sOutputFile:
    sOutputFile="Result_"+os.path.basename(sInputFile).split(".")[0]+".tsv"
    print("Warning : no outputfile -o defined, default {}".format(sOutputFile))

sExonListFile=options.exonlistfile
if not sExonListFile:
    sExonListFile="ExonList_"+os.path.basename(sInputFile).split(".")[0]+".txt"
    print("Warning : no exonlistfile -e defined, default {}".format(sExonListFile))

sMegablastResultFile=options.megablaststrand
bMegablast=True
if not sMegablastResultFile:
    bMegablast=False
    print("Warning : no megablaststrand -s defined, all Hsp keeped by default")

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
    
########################################################################
#MAIN
if __name__ == "__main__":
    if bMegablast:
        dRead2Strand=ExtractRead2Strand(sMegablastResultFile)
    
    oRefContent=TranscriptContent(sTargetId,sGffFile)
    print("Ref contains {} exons along {} isoformes".format(len(oRefContent.get_exonList()),len(oRefContent.get_transcript())))
    sExonList=oRefContent.describe_exon()
    WriteFile(sExonListFile,sExonList)
    
    if iTool==0:
        sTool="Blast"
        oFileContent=BlastContent(sInputFile)
    elif iTool==1:
        sTool="Yass"    
        oFileContent=YassContent(sInputFile)
        ##Yass Hack : sometimes, Yass can produce align already include into another align.
        ## the hack is made to remove this included alignment (interference into ExonCover and SelfCover computing)
        oYassCopyContent=copy.deepcopy(oFileContent)
        oFileContent.merge_withContent(oYassCopyContent)
        oFileContent.set_tool("yass")
    elif iTool==2:
        sTool="Exonerate"
        oFileContent=ExonerateContent(sInputFile,sRefFile,sFastaFile,dRead2Strand)
        #oFileContent.purge_alignList(["ch69_read19500_template_pass_BYK_CB_ONT_1_FAF04998_A"])
        #oFileContent.describe_AllHsp()
        oFileContent.check_divergence()
        #oFileContent.describe_AllHsp()
    print("{} length : {}nt".format(oFileContent.query_id,oFileContent.query_size))
    print("{} associated reads : {}".format(sTool,len(oFileContent.align_list)))
    ###TSV Ref2Read coord
    #WriteFile(sOutputFile,oFileContent.export_data())
    
    ##TSV exon covering file
    oFileContent.compute_self_cover_alignments()
    oFileContent.define_position_alignments()
    oFileContent.assign_exonCovering(oRefContent)
    sContent=oFileContent.describe_alignment_header(oRefContent,False)
    sContent+=oFileContent.describe_alignment(oRefContent,False)
    #print(sContent)
    WriteFile(sOutputFile,sContent)
    
    ###TSV Ref2Read coord
    WriteFile(".".join(sOutputFile.split(".")[:-1])+".paf",oFileContent.export_data())
