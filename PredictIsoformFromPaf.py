# coding: utf-8
"""Python3.6"""

import os
import sys
from optparse import OptionParser

from CurrentWorkflowClass import *

sCurrentVersionScript="v2"
########################################################################
'''
V2-2018/10/24
Use matrix logic

V1-2018/10/24
Load data from paf and predict exon and isoform

python PredictIsoformeFromPaf.py -p PAFFILE -t TARGETID -g GFFFILE 

PAFFILE: path to the paf file
TARGETID: gene of interest's EnsemblId
GFFFILE: path to the gff file of reference
'''
########################################################################
#Options
parser = OptionParser(conflict_handler="resolve")
parser.add_option("-p","--paffile", dest="paffile")
parser.add_option("-t","--targetid", dest="targetid")
parser.add_option("-g","--gfffile", dest="gfffile")
parser.add_option("-o","--outputfile", dest="outputfile")

(options, args) = parser.parse_args()

sInputFile=options.paffile
if not sInputFile:
    exit("Error : no paffile -p defined, process broken")

sOutputFile=options.outputfile
if not sOutputFile:
    sOutputFile="Result_"+os.path.basename(sInputFile).split(".")[0]+".tsv"
    print("Warning : no outputfile -o defined, default {}".format(sOutputFile))

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

########################################################################
#MAIN
if __name__ == "__main__":
    ##PREDICTION - Initialization
    oNanoporeContent=AlignedMatrixContent(sTargetId,sInputFile,sOutputFile)

    
    
