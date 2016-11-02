import sys
import operator

from sklearn.lda import LDA
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from sklearn import preprocessing


########################################################################
###Global variables
archFile        = ""
outputFile      = ""
colArchProtein  = 0
colArchDomain   = 1
groupFile       = ""
colGroupName    = 0
colGroupProtein = 1 
outputFile      = ""
sep             = "\t"

#Usage
usage = "python DA_score.py [options] -i inputFile -o outFile\nwhere basic options are:\n"
usage +="-v <x>\t: algorithm version, 1 for the simple version; 2 for changing centroids\n"

########################################################################
### Read parameters
def readParameters(args):
    global archFile, outputFile, colArchProtein, colArchDomain, groupFile, colGroupName, colGroupProtein; 
  
    for i in range(1,len(args)):
        if (args[i] == "-archFile"):
            archFile = args[i+1]
        elif (args[i] == "-colArchProtein"):
            colArchProtein = args[i+1]
            colArchProtein = int(colArchProtein)
        elif (args[i] == "-colArchDomain"):
            colArchDomain = args[i+1]
            colArchDomain =  int(colArchDomain)
        elif (args[i] == "-groupFile"):
            groupFile = args[i+1]
        elif (args[i] == "-colGroupName"):
            colGroupName = args[i+1]
            colGroupName =  int(colGroupName)
        elif (args[i] == "-colGroupProtein"):
            colGroupProtein = args[i+1]
            colGroupProtein =  int(colArchDomain)
        elif (args[i] == "-outputFile"):
            outputFile = args[i+1]
        elif (args[i] == "-h"):
            print (usage)
########################################################################
### Check parameters
def checkParameters():
    if (archFile == ""):
        print ("ERROR::Parameter -archFile is required\n")
        sys.exit(1);
    elif (outputFile == ""):
        print ("ERROR::Parameter -o outputFile is required\n")
        sys.exit(1);



########################################################################
###  readFileDict
def readFileDict(fileName, field, cols, sep):
    hashF = {}
    with open(fileName) as f:
        for line in f:
            line = line.strip()            
            arrayLine = line.split(sep) 
            key = arrayLine[field]          
            arrayLine = operator.itemgetter(cols)(arrayLine)            
            if (not hashF.has_key(key)):
                hashF[key] = []
            hashF[key].append(arrayLine)
                
    return hashF

########################################################################
###  readFileDict
def computeDAscore(hashGroup, hashArch):
    
    for group in hashGroup:
        listProt = hashGroup[group]
        for protein in listProt:
            if (hashArch.has_key(protein)):
                print hashArch[protein]
        print (group, listProt)
        sys.exit(1);
    
    
########################################################################
#################              MAIN                #####################
########################################################################

readParameters(sys.argv)
checkParameters()

hashGroup   = readFileDict(groupFile, colGroupName, colGroupProtein, sep);
hashArch    = readFileDict(archFile, colArchProtein, colArchDomain, sep);

computeDAscore(hashGroup, hashArch)
print (hashGroup)

