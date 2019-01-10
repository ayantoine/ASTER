# coding: utf-8
"""Python3.6"""

import os
import sys
from optparse import OptionParser
import random
import subprocess
import time
import datetime

sCurrentVersionScript="v11"
iTime1=time.time()
########################################################################
'''
V11-2018/12/20
Rework all the logic. Keep trace of starting alignment

python BlastWorkflow_Master.py -a ARGFILE

ARGFILE: txt file with all option
'''
########################################################################
#Options
parser = OptionParser(conflict_handler="resolve")
parser.add_option("-a","--argfile", dest="argfile")

(options, args) = parser.parse_args()

sArgFile=options.argfile
if not sArgFile:
    sys.exit("Error : no argfile -a defined, process broken")

########################################################################
#CONSTANT
#ROOTDIR=""
MASKDIR="wf_Maskdir"

READCOVER_THRESHOLD=""
GENECOVER_THRESHOLD=""

PYTHONDIR=""
SCRIPT1=""
SCRIPT2=""
SCRIPT3=""
SCRIPT4=""
SCRIPT5=""

GITDIR=""
SCRIPTGRAPH=""

BLASTDATABASE=""

ARGFILE_ASSIGNATION="="
ARGFILE_SEPARATOR=";"

AGREEDOWNLOAD=""

YEAR=datetime.date.today().year
MONTH=datetime.date.today().month
DAY=datetime.date.today().day

########################################################################
#CLASS
    
########################################################################
#Function 
def ParseArgFile(sFile):
    dListOfArgs={}
    for sLine in open(sFile):
        if ARGFILE_ASSIGNATION not in sLine:
            continue
        sLine=sLine.replace("\n","")
        tLine=sLine.split(ARGFILE_ASSIGNATION)
        sKey=tLine[0]
        sValue=tLine[1]
        if ARGFILE_SEPARATOR in sValue:
            tValue=sValue.split(ARGFILE_SEPARATOR)
            dListOfArgs[sKey]=tuple(tValue)
        else:
            dListOfArgs[sKey]=sValue
    return dListOfArgs

def LoadArgFile(sFile):
    dArgs=ParseArgFile(sFile)
    try:
        global READCOVER_THRESHOLD
        READCOVER_THRESHOLD=dArgs['READCOVER_THRESHOLD']
        global GENECOVER_THRESHOLD
        GENECOVER_THRESHOLD=dArgs['GENECOVER_THRESHOLD']
        
        global PYTHONDIR
        PYTHONDIR=dArgs['PYTHONDIR']
        global SCRIPT1
        SCRIPT1=PYTHONDIR+"/"+dArgs['SCRIPT1']
        global SCRIPT2
        SCRIPT2=PYTHONDIR+"/"+dArgs['SCRIPT2']
        global SCRIPT3
        SCRIPT3=PYTHONDIR+"/"+dArgs['SCRIPT3']
        global SCRIPT4
        SCRIPT4=PYTHONDIR+"/"+dArgs['SCRIPT4']
        global SCRIPT5
        SCRIPT5=PYTHONDIR+"/"+dArgs['SCRIPT5']
        
        global GITDIR
        GITDIR=dArgs['GITDIR']
        global SCRIPTGRAPH
        SCRIPTGRAPH=GITDIR+"/"+dArgs['SCRIPTGRAPH']
        
        global BLASTDATABASE
        BLASTDATABASE=dArgs['BLASTDATABASE']
        
        dbTarget=dArgs['GENEID']
    except KeyError as oError:
        sys.exit("Error : There is a missing value into the argsfile, process broken\n {}".format(oError))
    if isinstance(dbTarget,str):
        dbTarget=tuple([dbTarget])

    return dbTarget
 
def ExecuteBashCommand(scriptfile,bFirstLine=True):
    print("DEBUG:{}".format(scriptfile))
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
    
def CreateFolder(sFolder):
    if not os.path.isdir(sFolder):
        os.makedirs(sFolder)

########################################################################
#MAIN
if __name__ == "__main__":
    dbListOfTarget=LoadArgFile(sArgFile)
    tFolderList=[]
    
    for sGeneId in dbListOfTarget:
        try:
            print("WORKING ON {}".format(sGeneId))
            
            CreateFolder("./{}".format(sGeneId))
            tFolderList.append("./{}".format(sGeneId))
            
            print("--download...")
            ExecuteBashCommand("time python {1} -t {0} -g {0}/{0}.gff -f {0}/{0}.fa".format(
                sGeneId,SCRIPT1))
                    
            CreateFolder("{}/{}".format(sGeneId,MASKDIR))
            print("--repeatmasker...")
            ExecuteBashCommand("time RepeatMasker -pa 4 -s -species mouse -dir {0}/{1} {0}/{0}.fa".format(
                sGeneId,MASKDIR))
            if not os.path.isfile("{0}/{1}/{0}.fa.masked".format(sGeneId,MASKDIR)):
                print("WARNING 250 : No repetitive sequences present. Use .fa as .fa.masked")
                ExecuteBashCommand("cp {0}/{0}.fa {0}/{1}/{0}.fa.masked".format(
                    sGeneId,MASKDIR,FASTADIR))
            ExecuteBashCommand("mv {0}/{1}/{0}.fa.masked {0}/{0}.fa.masked".format(
                sGeneId,MASKDIR))
                
            ##Megablast pre-filter
            print("--megablast...")
            ExecuteBashCommand("time blastn -db {1} -query {0}/{0}.fa.masked -max_target_seqs 1000000 -outfmt 5 -out {0}/{0}.megablast.xml".format(
                sGeneId,BLASTDATABASE))
            print("--megablast parsing...")
            ExecuteBashCommand("time python {1} -b {0}/{0}.megablast.xml -t {0} -g {0}/{0}.gff -o {0}/{0}.megablast.tsv -e {0}/{0}.exon.txt".format(
                sGeneId,SCRIPT2))
            print("--extract read...")
            ExecuteBashCommand('time python {1} -r {3} -g {4} -f {2} -s {0}/{0}.megablast.tsv -o {0}/megablast.fa'.format(
                    sGeneId,SCRIPT3,BLASTDATABASE,READCOVER_THRESHOLD,GENECOVER_THRESHOLD))
            
            ##Exonerate
            print("--exonerate...")
            ExecuteBashCommand('time exonerate --showvulgar no --showtargetgff y --model est2genome {0}/megablast.fa {0}/{0}.fa.masked > {0}/{0}.exonerate.txt'.format(
                    sGeneId))
            print("--exonerate parsing...")
            ExecuteBashCommand('time python {1} -r {2} -b {0}/{0}.megablast.tsv -g {0}/{0}.fa.masked -e {0}/{0}.exonerate.txt -o {0}/{0}.exonerate.xml'.format(
                    sGeneId,SCRIPT4,BLASTDATABASE))
            print("--table...")
            ExecuteBashCommand("time python {1} -x {0}/{0}.exonerate.xml -r {0}/{0}.fa.masked -f {0}/megablast.fa -t {0} -g {0}/{0}.gff -o {0}/{0}.spliceSummary.tsv".format(
                sGeneId,SCRIPT5))
            
            
            #print("--exonerate parsing...")
            #ExecuteBashCommand("time python {1} -x {0}/{0}.exonerate.txt -t {0} -g {0}/{0}.gff -o {0}/{0}.exonerate.tsv -e {0}/{0}.exon.txt -f {0}/megablast.fa -r {0}/{0}.fa.masked -s {0}/{0}.megablast.tsv".format(
                #sGeneId,SCRIPT2))
                
            ###Graph
            #print("--graph...")
            #ExecuteBashCommand("time python {1} {0}/{0}.fa {0}/{0}.exonerate.paf {0}/{0}.exonerate.png".format(
                #sGeneId,SCRIPTGRAPH))
            #FILE=open("{0}/LaunchFile.sh".format(sGeneId),"w")
            #FILE.write("time python {1} ./{0}.fa ./{0}.exonerate.paf ./{0}.exonerate.png".format(
                #sGeneId,SCRIPTGRAPH))
            #FILE.close()
            
            ##Table
            #print("--table...")
            #ExecuteBashCommand("time python {1} -p {0}/{0}.exonerate.paf -r {0}/{0}.fa.masked -f {0}/megablast.fa -t {0} -g {0}/{0}.gff -o {0}/{0}.spliceSummary.tsv".format(
                #sGeneId,SCRIPT5))
            
            print("WORK ON {} DONE".format(sGeneId))
        except Exception as e:
            print("CRASH : {}".format(e))
        
    
    #print("GROUP FOLDERS")
    #sFolderPath="./{}analysisOn{}_{}{}{}".format(
                    #sCurrentVersionScript,len(dbListOfTarget),YEAR,MONTH,DAY)
    #CreateFolder(sFolderPath)
    #for sTargetFolder in tFolderList:
        #ExecuteBashCommand("mv ./{1} {0}/".format(sFolderPath,sTargetFolder))
    #print("DONE")

########################################################################    
iTime2=time.time()
iDeltaTime=iTime2-iTime1
print("Time : {}".format(iDeltaTime))
