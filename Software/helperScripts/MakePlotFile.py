'''
Created on Apr 7, 2014

@author: oliwa
'''

import argparse
import sys
from scriptutils import getFullPathOfURI, makeStringEndWith, makeStringNotEndWith
import glob as glob
import os as os

def createPlotString(postAnalyzerPath, description):
    """ Create the plotFile readable by the PostPlot*.py program. 
    
        Args:
            postAnalyzerPath: Path to the output from PostAnalyzer
            description: description of the results (for instance canonical or submatrix)
            
        Returns:
            plotFile readable by the PostPlot*.py program
    """
    analyzedFolders = glob.glob(postAnalyzerPath+"*/")
    scriptLine = []
    script = []
    for folder in analyzedFolders:
        if not (folder.endswith('_rigid/') or folder.endswith('_medium/') or folder.endswith('_difficult/') or folder.endswith('_all/')):
            continue 
        scriptLine = []
        scriptLine.append(folder)
        scriptLine.append(description)
        scriptLine.append("AssessmentWhole") # for now, later this needs to be refactored to made the argument optional
        difficulty = folder[folder.rfind("_")+1:-1]
        difficulty = difficulty.upper()[0]+difficulty[1:]
        scriptLine.append("Difficulty"+difficulty)
        script.append(scriptLine)
    # sort to have all, difficult, medium, rigid in this order (a before d before m before r)
    script.sort(key=lambda x: x[3], reverse=False)
    return script

def listOflistsToString(listOfLists):
    """ Convert a list of lists to a string, each list top level list separated by a newline.
    
        Args:
            listOfLists: a list of lists, containing strings in the lowest level
    
        Returns:
            listOfLists as a string, each top level list separated by a newline
     """
    output = []
    line = []
    for element in listOfLists:
        line = " ".join(element)
        line = line
        output.append(line)
    return "\n".join(output)

def main():
    """ Run the program MakePlotFile.py 
    
        Returns:
            path to the created plotFile.txt file
    """
    parser = argparse.ArgumentParser(description='This program creates a plot file for PostPlot*.py \n It assumes that the PostAnalysis has been already run (for example via the RunPostAnalyzer).',
                                     epilog='Result: The plot file is written to <resultsPath>_results/plotFile.txt')
    parser.add_argument('resultsPath', help='Path to the individual (per protein) results, to be plotted')
    parser.add_argument('optionalResultsPath', nargs='?', help='Optional path to the individual (per protein) canonical results, which can be used as comparison in the plots')
    parser.add_argument('--resultsDescriptor', help='A string describing the approach that created the results, default is \"submatrix\"')
    parser.add_argument('--optionalResultsDescriptor', help='A string describing the approach that created the optional results, default is \"canonical\"')
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()    
    
    # setup paths of individual and post analyzer results and the descriptors
    postAnalyzerSuffix = "_results/"
    
    if args.resultsDescriptor:
        resultsDescriptor = args.resultsDescriptor
    else:
        resultsDescriptor = "submatrix"
        
    if args.optionalResultsDescriptor:
        optionalResultsDescriptor = args.optionalResultsDescriptor
    else:
        optionalResultsDescriptor = "canonical"
    
    resultsPath = makeStringEndWith(getFullPathOfURI(args.resultsPath), "/")
    assert os.path.isdir(resultsPath)
    resultsPostAnalyzerPath = makeStringEndWith(makeStringNotEndWith(resultsPath, "/")+postAnalyzerSuffix, "/")
    assert os.path.isdir(resultsPostAnalyzerPath)
    
    plotString = createPlotString(resultsPostAnalyzerPath, resultsDescriptor)
    plotString = listOflistsToString(plotString)
    
    if args.optionalResultsPath:
        optionalResultsPath = makeStringEndWith(getFullPathOfURI(args.optionalResultsPath), "/")
        assert os.path.isdir(optionalResultsPath)
        optionalResultsPostAnalyzerPath = makeStringEndWith(makeStringNotEndWith(optionalResultsPath, "/")+postAnalyzerSuffix, "/")
        assert os.path.isdir(optionalResultsPostAnalyzerPath)
        
        optionalPlotString = createPlotString(optionalResultsPostAnalyzerPath, optionalResultsDescriptor)
        optionalPlotString = listOflistsToString(optionalPlotString)
        plotString = plotString + "\n" + optionalPlotString

    with open(resultsPostAnalyzerPath+"plotFile.txt", "w") as text_file:
        text_file.write(plotString)

    print getFullPathOfURI(resultsPostAnalyzerPath+"plotFile.txt")
    return getFullPathOfURI(resultsPostAnalyzerPath+"plotFile.txt")
if __name__ == '__main__':
    main()