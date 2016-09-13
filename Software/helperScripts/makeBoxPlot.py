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

def returnSubsetOflist(list1, indexToAppend):
    subsetList = []
    for idx in indexToAppend:
        subsetList.append(list1[idx])
    return np.array(subsetList)

def relativeReductions(cNMA, NMA):
    return ( (cNMA/NMA) - 1.0) * 100

def main():
    parser = argparse.ArgumentParser(description='Make a box plot of values, such as RMSD Reduction values')
    parser.add_argument('resultsPath', help='Folder with subfolders, each having protein results')
    parser.add_argument('reductionFile', help='Name of the file containing the reduction values')
    parser.add_argument('originalValue', help='Initial value, that the experiments where reducing from')
    parser.add_argument('outputFile', help='Output file')
    parser.add_argument('-category', help='Only use proteins that foldernames match the names in this list')
    parser.add_argument('-title', help='Add custom title to the plot')
    parser.add_argument('-relativeReduction', help='Plot relative reduction compared against this folder with protein results')
    
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
        
    # reductionsListOfLists
    reductionsListOfLists = []
    
    listOfproteinPaths = sorted(glob.glob(args.resultsPath+"*/"))
   
    # plot variables
    fontSize = 19
    
    if args.title:
        plotTitle = args.title
    else:
        plotTitle = ""

    for proteinPath in listOfproteinPaths:
        # skip if protein in not in category
        if args.category:
            pdbName = np.loadtxt(proteinPath+"pdbName.txt", dtype="string")
            pdbName = np.array_str(pdbName)
            if pdbName not in categoryProteins:
                continue
        #print proteinPath+args.reductionFile
        assert os.path.isfile(proteinPath+args.reductionFile)
        assert os.path.isfile(proteinPath+args.originalValue)
        
        reductions = np.loadtxt(proteinPath+args.reductionFile)
        originalValue = np.loadtxt(proteinPath+args.originalValue)
        
        # make RMSD reduction list
        reductionsPercentage = makeReductionList(originalValue, reductions)
        
        # make subset of only certain x axis steppoints
        indexToAppend = [9, 19, 20, 21, 22, 23, 24, 25]
        reductionsPercentageSubset = returnSubsetOflist(reductionsPercentage, indexToAppend)
        
        # create relative reductions
        
        if args.relativeReduction:
            assert os.path.isfile(args.relativeReduction+pdbName+"/"+args.reductionFile)
            assert os.path.isfile(args.relativeReduction+pdbName+"/"+args.originalValue)
            reductions_toCompareWith = np.loadtxt(args.relativeReduction+pdbName+"/"+args.reductionFile)
            originalValue_toCompareWith = np.loadtxt(args.relativeReduction+pdbName+"/"+args.originalValue)
            assert np.isclose(originalValue, originalValue_toCompareWith)           
            
            # make RMSD reduction list
            reductionsPercentage_toCompareWith = makeReductionList(originalValue_toCompareWith, reductions_toCompareWith)
            reductionsPercentageSubset_toCompareWith = returnSubsetOflist(reductionsPercentage_toCompareWith, indexToAppend)
            reductionsPercentageSubset = relativeReductions(reductionsPercentageSubset, reductionsPercentageSubset_toCompareWith)
        
        # append to the overall reductionListOfLists
        reductionsListOfLists.append(reductionsPercentageSubset)

    reductionsListOfLists = np.array(reductionsListOfLists)
    
    # box plot 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    reductionsListOfLists_T = reductionsListOfLists.T
    print "reductionsListOfLists.shape: ", reductionsListOfLists.shape
    
    xAxis = ['10', '20', '30', '40', '50', '60', '70', '80']
    
    ax.set_xticklabels(xAxis)
    ax.set_title("RMSD red. % "+plotTitle, size=fontSize-4)
    ax.set_xlabel("Modes", size=fontSize)
    ax.set_ylabel("RMSD Reduction (%)", size=fontSize) 
#     ax.boxplot([
#                 reductionsListOfLists_T[0],
#                 reductionsListOfLists_T[1],
#                 reductionsListOfLists_T[2],
#                 reductionsListOfLists_T[3],
#                 reductionsListOfLists_T[4],
#                 reductionsListOfLists_T[5],
#                 reductionsListOfLists_T[6],
#                 reductionsListOfLists_T[7],
#                 ])
    
    ax.plot([10, 20, 30, 40, 50, 60, 70, 80], reductionsListOfLists_T.mean(axis=1))
    
    plt.legend(loc='best', prop={'size':12})
    plt.savefig(args.outputFile+'.eps', bbox_inches='tight', dpi=1200)
    plt.savefig(args.outputFile+'.pdf', bbox_inches='tight', dpi=1200)
    
if __name__ == '__main__':
    main()