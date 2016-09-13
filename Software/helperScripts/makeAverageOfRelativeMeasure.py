'''
Created on Sep 26, 2014

@author: oliwa
'''

import argparse
import sys
import os
import numpy as np
from prody import *
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def percentageDecrease(x, y):
    """Return the percentage decrease from x towards y
    
        Args:
            x: original value
            y: changed value
        
        Returns: percentage decrease from x towards y.
            
    """
    return (float(x)-float(y))/float(x)

def makeReductionList(originalValue, absoluteReductions):
    reductionsPercentage = []
    for element in absoluteReductions:
        reductionsPercentage.append(percentageDecrease(originalValue, element))
    return (np.array(reductionsPercentage) * 100)

def skipIfNAN(proteinPath):
    """ Test if there is a NAN (not a number) in the lists """
    overlapArrayWhole = None
    overlapArrayInterface = None
    overlapTApproxWhole = None
    overlapTApproxInterface = None
    try:
        overlapArrayWhole = np.loadtxt(proteinPath+"overlapArrayWhole.txt")
    except IOError:
        pass
    try:
        overlapArrayInterface = np.loadtxt(proteinPath+"overlapArrayInterface.txt")
    except IOError:
        pass            
    try:
        overlapTApproxWhole = np.loadtxt(proteinPath+"overlapTApproxWhole.txt")
    except IOError:
        pass            
    try:
        overlapTApproxInterface = np.loadtxt(proteinPath+"overlapTApproxInterface.txt")
    except IOError:
        pass            
    if overlapArrayWhole is not None and np.isnan(overlapArrayWhole).any():
        #print "skipped"
        return True
    if overlapArrayInterface is not None and np.isnan(overlapArrayInterface).any():
        #print "skipped"
        return True
    if overlapTApproxWhole is not None and np.isnan(overlapTApproxWhole).any():
        #print "skipped"
        return True
    if overlapTApproxInterface is not None and np.isnan(overlapTApproxInterface).any():
        #print "skipped"
        return True
    return False

def returnSubsetOflist(list1, indexToAppend):
    subsetList = []
    for idx in indexToAppend:
        subsetList.append(list1[idx])
    return np.array(subsetList)

def relativeReductions(cNMA, NMA):
    return ( (cNMA/NMA) - 1.0) * 100

def relativeReductionsExclude(cNMA, NMA, excludeBelowPercentage):
    assert len(cNMA) == len(NMA)
    result = []
    for cNMA_value, NMA_value in zip(cNMA, NMA):
        if (cNMA_value < excludeBelowPercentage) and (NMA_value < excludeBelowPercentage):
            relativeImprovement = 0.0
        else:
            relativeImprovement = ( (cNMA_value/NMA_value) - 1.0) * 100
        result.append(relativeImprovement)
    return np.array(result)

def main():
    parser = argparse.ArgumentParser(description='Calculate and output relative measures of RMSD Reductions between two methods')
    parser.add_argument('resultsPath', help='Folder with subfolders, each having protein results')
    parser.add_argument('reductionFile', help='Name of the file containing the reduction values')
    parser.add_argument('originalValue', help='Initial value, that the experiments where reducing from')
    parser.add_argument('outputFile', help='Output file, filename excluding extension')
    parser.add_argument('-outputPath', help='Output path')
    parser.add_argument('-category', help='Only use proteins that are in foldernames which match the names in this list')
    parser.add_argument('-relativeReduction', help='Plot relative reduction compared against this folder with protein results')
    parser.add_argument('-excludeBelowXPercent', help='If both reductions are below the percentage X, set relative improvement to 0.')
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()    

    assert os.path.isdir(args.resultsPath)
    if args.category:
        assert os.path.isfile(args.category)
        categoryProteins = set(list((np.loadtxt(args.category, dtype="string"))))
        
    if args.relativeReduction:
        assert os.path.isdir(args.relativeReduction)
        
    if args.excludeBelowXPercent:
        excludeBelowPercentage = float(args.excludeBelowXPercent)
        print "excluding below ", excludeBelowPercentage," percent for both percentage improvements"

    # reductionsListOfLists
    reductionsListOfLists = []
    
    # tocomparewithList
    reductionsListOfLists_toCompareWith = []
    
    # relativeMeasure_perProtein
    relativeMeasure_perProtein = []
    
    listOfproteinPaths = sorted(glob.glob(args.resultsPath+"*/"))
    
    counter = 0
    for proteinPath in listOfproteinPaths:
        pdbName = os.path.basename(os.path.normpath(proteinPath))
        # skip if protein in not in category
        if args.category:
            if os.path.isfile(proteinPath+"pdbName.txt"):
                assert pdbName == np.array_str(np.loadtxt(proteinPath+"pdbName.txt", dtype="string"))
            if pdbName not in categoryProteins:
                continue
        # skip if NaNs present, to test
        if skipIfNAN(proteinPath):
            continue
        proteinInterfacePath = proteinPath.replace("_whole_", "_interface_")
        if skipIfNAN(proteinInterfacePath):
            continue            

        # create relative reductions
        if args.relativeReduction:
            assert os.path.isfile(args.relativeReduction+pdbName+"/"+args.reductionFile)
            assert os.path.isfile(args.relativeReduction+pdbName+"/"+args.originalValue)
            reductions_toCompareWith = np.loadtxt(args.relativeReduction+pdbName+"/"+args.reductionFile)
            originalValue_toCompareWith = np.loadtxt(args.relativeReduction+pdbName+"/"+args.originalValue)
            assert np.isclose(np.loadtxt(proteinPath+args.originalValue), originalValue_toCompareWith)           
             
            # make RMSD reduction list
            reductionsPercentage_toCompareWith = makeReductionList(originalValue_toCompareWith, reductions_toCompareWith)
            reductionsPercentage_toCompareWith = reductionsPercentage_toCompareWith[:26]
            reductionsListOfLists_toCompareWith.append(reductionsPercentage_toCompareWith)

        #print proteinPath+args.reductionFile
        assert os.path.isfile(proteinPath+args.reductionFile)
        assert os.path.isfile(proteinPath+args.originalValue)
        reductions = np.loadtxt(proteinPath+args.reductionFile)
        originalValue = np.loadtxt(proteinPath+args.originalValue)
        
        # make RMSD reduction list
        reductionsPercentage = makeReductionList(originalValue, reductions)
        reductionsPercentage = reductionsPercentage[:26]
        reductionsListOfLists.append(reductionsPercentage)
        
        if args.relativeReduction:
            if not args.excludeBelowXPercent:
                reductionsPercentageSubset = relativeReductions(reductionsPercentage, reductionsPercentage_toCompareWith)
            else:
                reductionsPercentageSubset = relativeReductionsExclude(reductionsPercentage, reductionsPercentage_toCompareWith, excludeBelowPercentage)
#             ###
#             print pdbName
#             print "original RMSD: ", originalValue
#             print "reductions cNMA  : ", reductions[:26].round(3)
#             print "reductions cNMA %: ", reductionsPercentage.round(3) 
#             print "reductions NMA   : ", reductions_toCompareWith[:26].round(3)
#             print "reductions NMA % : ", reductionsPercentage_toCompareWith.round(3)
#             print "relative measure : ", reductionsPercentageSubset.round(3)
#             ###
            relativeMeasure_perProtein.append(reductionsPercentageSubset)

        counter = counter + 1
    # make the mean array
    reductionsListOfLists = np.array(reductionsListOfLists)
    reductionMeans = reductionsListOfLists.mean(axis=0)
    
    if args.relativeReduction:
        # make the mean array
        reductionsListOfLists_toCompareWith = np.array(reductionsListOfLists_toCompareWith)
        reductionsListOfLists_toCompareWith_Means = reductionsListOfLists_toCompareWith.mean(axis=0)        
        relativeMeasure = relativeReductions(reductionMeans, reductionsListOfLists_toCompareWith_Means)
        #np.savetxt("relativeMeasure"+args.outputFile+".txt", relativeMeasure, fmt="%s")
        
        relativeMeasure_perProtein = np.array(relativeMeasure_perProtein)
        relativeMeasure_perProtein_Means = relativeMeasure_perProtein.mean(axis=0)
        
        if args.outputPath:
            assert os.path.isdir(args.outputPath)
            np.savetxt(args.outputPath+args.outputFile+".txt", relativeMeasure_perProtein_Means, fmt="%s", header=args.outputFile+"_"+str(counter))
        else:
            np.savetxt(args.outputFile+".txt", relativeMeasure_perProtein_Means, fmt="%s", header=args.outputFile+"_"+str(counter))
        #np.savetxt("toCompareWith_Means"+args.outputFile+".txt", reductionsListOfLists_toCompareWith_Means, fmt="%s")

    # output the mean
    #np.savetxt("reductionMeans_"+args.outputFile+".txt", reductionMeans, fmt="%s")

if __name__ == '__main__':
    main()