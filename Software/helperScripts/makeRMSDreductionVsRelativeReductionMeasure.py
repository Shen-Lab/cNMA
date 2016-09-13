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
        print "skipped"
        return True
    if overlapArrayInterface is not None and np.isnan(overlapArrayInterface).any():
        print "skipped"
        return True
    if overlapTApproxWhole is not None and np.isnan(overlapTApproxWhole).any():
        print "skipped"
        return True
    if overlapTApproxInterface is not None and np.isnan(overlapTApproxInterface).any():
        print "skipped"
        return True
    return False

def returnSubsetOflist(list1, indexToAppend):
    subsetList = []
    for idx in indexToAppend:
        subsetList.append(list1[idx])
    return np.array(subsetList)

def relativeReductions(cNMA, NMA):
    return ( (cNMA/NMA) - 1.0) * 100

def main():
    parser = argparse.ArgumentParser(description='Make a plot with the relative RMSD reduction on the X-axis, and improvement (relative reduction) on the Y-axis. The output are columnwise related files with xaxis, yaxis, pdbname and class being in the same column.')
    parser.add_argument('resultsPath', help='Folder with subfolders, each having protein results')
    parser.add_argument('reductionFile', help='Name of the file containing the reduction values')
    parser.add_argument('originalValue', help='Initial value, that the experiments where reducing from')
    parser.add_argument('categoryPath', help='Path to the categories in txt files all.txt difficult.txt medium.txt rigid.txt')
    parser.add_argument('relativeReduction', help='Plot relative reduction compared against this folder with protein results. For instance, this folder could have the conventional NMA results.')
    parser.add_argument('k', help='k, at which array index position of the stepsize of the expanding superset (number of modes used for reduction) to take the results')
    parser.add_argument('-outputPath', help='Output path')
    parser.add_argument('-header', help='header of output files')
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()    

    assert os.path.isdir(args.resultsPath)
    assert os.path.isdir(args.categoryPath)
    categoryProteins_difficult = set(list((np.loadtxt(args.categoryPath+"difficult.txt", dtype="string"))))
    categoryProteins_medium = set(list((np.loadtxt(args.categoryPath+"medium.txt", dtype="string"))))
    categoryProteins_rigid = set(list((np.loadtxt(args.categoryPath+"rigid.txt", dtype="string"))))
    assert os.path.isdir(args.relativeReduction)
        
    if args.header:
        header = str(args.header)
    else:
        header = ""
        
    k = int(args.k)
        
    # reductionsListOfLists
    reductionsListOfLists = []
    
    # tocomparewithList
    reductionsListOfLists_toCompareWith = []
    
    # relativeMeasure_perProtein
    relativeMeasure_perProtein = []
    
    listOfproteinPaths = sorted(glob.glob(args.resultsPath+"*/"))
    
    counter = 0
    
    #variables to output
    xAxis = []
    yAxis = []
    pdbNames = []
    categories = []
    
    for proteinPath in listOfproteinPaths:
        currentCategory = None
        
        pdbName = os.path.basename(os.path.normpath(proteinPath))
        # skip if protein in not in category
        if args.categoryPath:
            if os.path.isfile(proteinPath+"pdbName.txt"):
                assert pdbName == np.array_str(np.loadtxt(proteinPath+"pdbName.txt", dtype="string"))
            if pdbName in categoryProteins_difficult:
                currentCategory = "difficult"
            elif pdbName in categoryProteins_medium:
                currentCategory = "medium"
            elif pdbName in categoryProteins_rigid:
                currentCategory = "rigid"
            else:
                print "Protein", pdbName ," not found in category"
                raise Exception("Protein not found in category")
            
        # skip if NaNs present, to test
        if skipIfNAN(proteinPath):
            continue
        proteinInterfacePath = proteinPath.replace("_whole_", "_interface_")
        if skipIfNAN(proteinInterfacePath):
            continue            

        # create relative reductions
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
        xAxis.append(reductionsPercentage[k])
        reductionsPercentage = reductionsPercentage[:26]
        reductionsListOfLists.append(reductionsPercentage)

        reductionsPercentageSubset = relativeReductions(reductionsPercentage, reductionsPercentage_toCompareWith)
        yAxis.append(reductionsPercentageSubset[k])
        pdbNames.append(pdbName)
        categories.append(currentCategory)
#         ###
#         print pdbName
#         print "original RMSD: ", originalValue
#         print "reductions cNMA  : ", reductions
#         print "reductions cNMA %: ", reductionsPercentage 
#         print "reductions NMA   : ", reductions_toCompareWith
#         print "reductions NMA % : ", reductionsPercentage_toCompareWith
#         print "relative measure : ", reductionsPercentageSubset
#         ###
        relativeMeasure_perProtein.append(reductionsPercentageSubset)

        counter = counter + 1
        print pdbName
        
    if args.outputPath:
        assert os.path.isdir(args.outputPath)
        np.savetxt(args.outputPath+"xAxis_"+str(k)+".txt", xAxis, fmt="%.3f", header=header+"_xaxis_modeIndex_"+str(k)+"_"+str(counter))
        np.savetxt(args.outputPath+"yAxis_"+str(k)+".txt", yAxis, fmt="%.3f", header=header+"_yaxis_modeIndex_"+str(k)+"_"+str(counter))
        np.savetxt(args.outputPath+"pdbNames_"+str(k)+".txt", pdbNames, fmt="%s", header=header+"_pdbName_"+str(counter))
        np.savetxt(args.outputPath+"categories_"+str(k)+".txt", categories, fmt="%s", header=header+"_categories_"+str(counter))
    else:
        np.savetxt("xAxis_"+str(k)+".txt", xAxis, fmt="%.3f", header=header+"_xaxis_modeIndex_"+str(k)+"_"+str(counter))
        np.savetxt("yAxis_"+str(k)+".txt", yAxis, fmt="%.3f", header=header+"_yaxis_modeIndex_"+str(k)+"_"+str(counter))
        np.savetxt("pdbNames_"+str(k)+".txt", pdbNames, fmt="%s", header=header+"_pdbName_"+str(counter))
        np.savetxt("categories_"+str(k)+".txt", categories, fmt="%s", header=header+"_categories_"+str(counter))

if __name__ == '__main__':
    main()