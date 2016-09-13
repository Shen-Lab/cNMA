'''
Created on Apr 2, 2014

@author: oliwa
'''

import os as os
import sys as sys
import traceback as traceback
from scriptutils import mkdir_p, makeStringEndWith, makeStringNotEndWith
import argparse as argparse 

def main():
    """ Run the RunPostAnalyzer.py program """
    parser = argparse.ArgumentParser(description="Script to run the PostAnalyzer on a folder for the difficulty categorizations: all, difficult, medium and rigid and output the results a results folder. The absolute path to the PostAnalysis.py program can be specified in the variable \"pathToPostAnalysis\"", 
                                     epilog="Returns: a folder <whole>_results/ with the postAnalyzer results.\n Example: $ python RunPostAnalyzer.py --onlypdb --whole /home/oliwa/nmaData/output/testing/complexresults/ --interface /home/oliwa/nmaData/output/testing/complexresultsinterface/ --subset /home/oliwa/nmaData/output/testing/mysubset.txt -> the output is the folder <../whole_results> with individual assessment folder results inside, one per difficulty.")
    parser.add_argument('--whole', help='Absolute folderpath with whole results per protein')
    parser.add_argument('--interface', help='Absolute folderpath with interface results per protein')
    parser.add_argument('--subset', help='Absolute folderpath with benchmark difficulty categorizations in <subset>/all.txt <subset>/difficult.txt <subset>/medium.txt <subset>/rigid.txt')
    parser.add_argument('--postAnalysisProgramPath', help='Path to the PostAnalysis program, default is \"/home/oliwa/workspace/TNMA1/src/PostAnalysis.py\"')
    parser.add_argument('--pythonPath', help='Use Absolute path to the python interpreter, for example /home/oliwa/anaconda/bin/python')
    parser.add_argument('--onlypdb', action="store_true", help="Only consider the first 4 letters of the foldernames when taking a subset")
    parser.add_argument('-excludeAfterDistance', help='exclude proteins that have their x-measure (distance) bigger than this value')    
    parser.add_argument('--excludeBasedonDistancePath', help="Exclude proteins from analysis if their distances are above a threshold.")    
    parser.add_argument('--excludefirstKModes', help="Exclude first k modes from overlap, correlation and collectivity analysis (for example if the first 6 modes are trivial normal modes, give the value 6).")    
    parser.add_argument('--excludeFromList', help="Path to list with protein names to only use (exclude the others) from postanalysis. ")        
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    
    if args.pythonPath:
        pythonPath = args.pythonPath
    else:
        pythonPath = "python"    
        
    
    #inputString to the PostAnalysis
    inputString = ""
    
    if args.whole:
        assert os.path.isdir(args.whole)
        postAnalysisFolderPath = makeStringEndWith(makeStringNotEndWith(args.whole, "/")+"_results/", "/")
        inputResultsPathWithoutSlash = makeStringNotEndWith(args.whole, "/")
        inputString = "--whole "+args.whole
    if args.interface:
        assert os.path.isdir(args.interface)
        if not args.whole:
            postAnalysisFolderPath = makeStringEndWith(makeStringNotEndWith(args.interface, "/")+"_results/", "/")
            inputResultsPathWithoutSlash = makeStringNotEndWith(args.interface, "/")
        inputString += " --interface "+args.interface
    assert os.path.isdir(args.subset)
    
    if args.postAnalysisProgramPath:
        assert os.path.isfile(args.postAnalysisProgramPath)
        pathToPostAnalysis = args.postAnalysisProgramPath
    else:
        pathToPostAnalysis = "/home/oliwa/workspace/TNMA1/src/PostAnalysis.py"
    
    benchmarkPath = makeStringEndWith(args.subset, "/")
    
    mkdir_p(postAnalysisFolderPath)
    assert os.path.isdir(postAnalysisFolderPath)
    
    # path variables to benchmark difficulties
    allPathToFile = benchmarkPath+"all.txt"
    difficultPathToFile = benchmarkPath+"difficult.txt"
    mediumPathToFile = benchmarkPath+"medium.txt"
    rigidPathToFile = benchmarkPath+"rigid.txt"
    
    assert os.path.isfile(difficultPathToFile)
    assert os.path.isfile(mediumPathToFile)
    assert os.path.isfile(rigidPathToFile)
    
    # Run PostAnalysis.py directly via shellcommands 
    if args.onlypdb:
        try:
            shellCommand = pythonPath+" "+pathToPostAnalysis+" --onlypdb "+inputString+ " --subset " + allPathToFile
            os.system(shellCommand)
            shellCommand = pythonPath+" "+pathToPostAnalysis+" --onlypdb "+inputString + " --subset " + difficultPathToFile
            os.system(shellCommand)
            shellCommand = pythonPath+" "+pathToPostAnalysis+" --onlypdb "+inputString + " --subset " + mediumPathToFile
            os.system(shellCommand)
            shellCommand = pythonPath+" "+pathToPostAnalysis+" --onlypdb "+inputString + " --subset " + rigidPathToFile
            os.system(shellCommand)
            shellCommand = "mv "+inputResultsPathWithoutSlash+"_assessment* "+postAnalysisFolderPath
            os.system(shellCommand)
        except Exception, err:
            print "Exception occurred when running the RunPostAnalyzer with os.system(shellCommand): ", err
            print traceback.format_exc()
    else:
        try:
            shellCommand = pythonPath+" "+pathToPostAnalysis+" "+inputString+ " --subset " + allPathToFile
            if args.excludeBasedonDistancePath:
                shellCommand = shellCommand + " -excludeAfterDistance " + args.excludeAfterDistance + " --excludeBasedonDistancePath " + args.excludeBasedonDistancePath
            if args.excludefirstKModes:
                shellCommand = shellCommand + " --excludefirstKModes " + args.excludefirstKModes
            if args.excludeFromList:    
                shellCommand = shellCommand + " --excludeFromList " + args.excludeFromList
            os.system(shellCommand)
            
            shellCommand = pythonPath+" "+pathToPostAnalysis+" "+inputString + " --subset " + difficultPathToFile
            if args.excludeBasedonDistancePath:
                shellCommand = shellCommand + " -excludeAfterDistance " + args.excludeAfterDistance + " --excludeBasedonDistancePath " + args.excludeBasedonDistancePath            
            if args.excludefirstKModes:
                shellCommand = shellCommand + " --excludefirstKModes " + args.excludefirstKModes            
            if args.excludeFromList:    
                shellCommand = shellCommand + " --excludeFromList " + args.excludeFromList            
            os.system(shellCommand)
            
            shellCommand = pythonPath+" "+pathToPostAnalysis+" "+inputString + " --subset " + mediumPathToFile
            if args.excludeBasedonDistancePath:
                shellCommand = shellCommand + " -excludeAfterDistance " + args.excludeAfterDistance + " --excludeBasedonDistancePath " + args.excludeBasedonDistancePath            
            if args.excludefirstKModes:
                shellCommand = shellCommand + " --excludefirstKModes " + args.excludefirstKModes            
            if args.excludeFromList:    
                shellCommand = shellCommand + " --excludeFromList " + args.excludeFromList            
            os.system(shellCommand)
            
            shellCommand = pythonPath+" "+pathToPostAnalysis+" "+inputString + " --subset " + rigidPathToFile
            if args.excludeBasedonDistancePath:
                shellCommand = shellCommand + " -excludeAfterDistance " + args.excludeAfterDistance + " --excludeBasedonDistancePath " + args.excludeBasedonDistancePath            
            if args.excludefirstKModes:
                shellCommand = shellCommand + " --excludefirstKModes " + args.excludefirstKModes            
            if args.excludeFromList:    
                shellCommand = shellCommand + " --excludeFromList " + args.excludeFromList            
            os.system(shellCommand)
            
            shellCommand = "mv "+inputResultsPathWithoutSlash+"_assessment* "+postAnalysisFolderPath
            os.system(shellCommand)
        except Exception, err:
            print "Exception occurred when running the RunPostAnalyzer with os.system(shellCommand): ", err
            print traceback.format_exc()
    print postAnalysisFolderPath
    return postAnalysisFolderPath
if __name__ == '__main__':
    main()