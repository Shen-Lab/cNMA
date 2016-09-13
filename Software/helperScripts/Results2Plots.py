'''
Created on Apr 17, 2014

@author: oliwa
'''

"""
This program runs several helper programs to go from individual NMA results to
analysis and plots
"""


from scriptutils import runBashCommand, getAllAfterLastChar, makeStringEndWith, getFullPathOfURI
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='This program runs several helper programs to go from individual NMA results to analysis and plots',
                                     epilog='Result: The results are written to <uniquePrefix_Path>_results/')
    parser.add_argument('--uniqueResultsPrefix_Path', help='Absolute path and unique prefix of the individual result file folders, for example /home/oliwa/nmaData/output/2bb_Individual_U1_results/test2/2bbb_Individual_k05_A10_lambdaC_enforceT_U1')
    parser.add_argument('--uniqueResultsPrefixInterface_Path', help='Absolute path and unique prefix of the individual interface result file folders, for example /home/oliwa/nmaData/output/2bb_Individual_U1_results/test2/2bbb_Individual_k05_A10_lambdaC_enforceT_U1')
    parser.add_argument('--resultsDescriptor', help='Short string descriptor of method used, this will be used in plot legends, for example "U1_submatrix"')
    parser.add_argument('--canonicalResults_Path', nargs='?', help='Absolute path to the individual canonical results, for example /home/oliwa/nmaData/output/2bb_Individual_U1_results/2bbb_Individual_canonical/"')    
    parser.add_argument('--canonicalResultsInterface_Path', nargs='?', help='Absolute path to the individual interface canonical results, for example /home/oliwa/nmaData/output/2bb_Individual_U1_results/2bbb_Individual_canonical/"')    
    # optional flags
    parser.add_argument('--categoryMultiplier', help='Integer that is a multiplier of the category size for success failure plotrates. For 2c, set this to 10, default is 1.')
    parser.add_argument('--dontCombine', action="store_true", help='Take the first argument as the individual results path.')
    parser.add_argument('--dontCompare', action="store_true", help='Do not compare results with optional (for instance canonical) results, just analyze and plot them.')
    parser.add_argument('--isComplex', action="store_true", help='Results are from complex NMA data, not individual protein, set benchmark categories accordingly.')
    parser.add_argument('--pythonPath', help='Use Absolute path to the python interpreter, for example /home/oliwa/anaconda/bin/python')
    parser.add_argument('-excludeAfterDistance', help='exclude proteins that have their x-measure (distance) bigger than this value')
    parser.add_argument('--excludeBasedonDistancePath', help="Exclude proteins from analysis if their distances are above a threshold.")
    parser.add_argument('--sf_benchmarkSubsetPath', help="path to the benchmark difficulty categories for the sf rates")
    parser.add_argument('--excludefirstKModes', help="Exclude first k modes from overlap, correlation and collectivity analysis (for example if the first 6 modes are trivial normal modes, give the value 6).")    
    parser.add_argument('--plotOnlyReduction', action="store_true", help='If set, only RMSD reduction and associated plots are made, not measures plots/barplots (overlap, collectivity etc.)')
    parser.add_argument('--excludeFromList', help="Path to list with protein names to only use (exclude the others) from postanalysis. ")    
    
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    
    if args.categoryMultiplier:
        categoryMultiplier = args.categoryMultiplier
    else:
        categoryMultiplier = 1.0    
      
    # exclude first k modes if argument has been given
    if args.excludefirstKModes:
        excludefirstKModes = int(args.excludefirstKModes)
    else:
        excludefirstKModes = None        
    
    #paths
    if args.isComplex:                                                                                      # /allinterfaceSuperposedPdbs/ or /benchmark40/ or /eigenspectrumCases/
        benchmarkSubset_Path="/home/oliwa/workspace/TNMA1/src/BenchmarkAssessmentsOfDifficultyBoundComplex/benchmark40/"
        sf_benchmarkSubsetPath = benchmarkSubset_Path
        taxonomyPath = "/home/oliwa/workspace/TNMA1/src/BenchmarkAssessmentsOfDifficulty/enzymeantibody/complexonlypdbs/"
    else:
        benchmarkSubset_Path="/home/oliwa/workspace/TNMA1/src/BenchmarkAssessmentsOfDifficulty/allinterfaceSuperposed/"
        sf_benchmarkSubsetPath = benchmarkSubset_Path
        taxonomyPath = "/home/oliwa/workspace/TNMA1/src/BenchmarkAssessmentsOfDifficulty/enzymeantibody/with_r_u_suffix/"
    
    if args.sf_benchmarkSubsetPath:
        sf_benchmarkSubsetPath = args.sf_benchmarkSubsetPath
        
    if args.pythonPath:
        pythonPath = args.pythonPath
    else:
        pythonPath = "python"
    AnalysisPlotter_Path="/home/oliwa/workspace/TNMA1/src/AnalysisPlotter.py"
    SuccessFailureAnalyzer_Path = "/home/oliwa/workspace/TNMA1/src/SuccessFailureAnalyzer.py"
    SuccessFailurePlotter_Path = "/home/oliwa/workspace/TNMA1/src/SuccessFailurePlotter.py"
    VarargsPlotter_Path = "/home/oliwa/workspace/TNMA1/src/VarargsPlotter.py"
    
    print "Benchmark difficulties in:", benchmarkSubset_Path
    print "Sf rates difficulties  in:", sf_benchmarkSubsetPath
    
    runPostAnalyzerInput = ""
    if not args.dontCombine:
        # ResultsCombiner.py
        if args.uniqueResultsPrefix_Path:
            bashCommand = pythonPath+" ResultsCombiner.py "+args.uniqueResultsPrefix_Path
            combinedIndividualResults_Path = runBashCommand(bashCommand).strip()
            runPostAnalyzerInput += "--whole "+combinedIndividualResults_Path
        if args.uniqueResultsPrefixInterface_Path:
            bashCommand = pythonPath+" ResultsCombiner.py "+args.uniqueResultsPrefixInterface_Path
            combinedIndividualResultsInterface_Path = runBashCommand(bashCommand).strip()
            runPostAnalyzerInput += " --interface "+combinedIndividualResultsInterface_Path        
    else:
        # Take the first argument as the combined results path
        if args.uniqueResultsPrefix_Path:
            combinedIndividualResults_Path = makeStringEndWith(getFullPathOfURI(args.uniqueResultsPrefix_Path), "/")
            runPostAnalyzerInput += "--whole "+combinedIndividualResults_Path
        if args.uniqueResultsPrefixInterface_Path:
            combinedIndividualResultsInterface_Path = makeStringEndWith(getFullPathOfURI(args.uniqueResultsPrefixInterface_Path), "/")
            runPostAnalyzerInput += " --interface "+combinedIndividualResultsInterface_Path        

    print combinedIndividualResults_Path

    # RunPostAnalyzer.py
    if args.isComplex:
        bashCommand = pythonPath+" RunPostAnalyzer.py --pythonPath "+pythonPath+" --onlypdb "+runPostAnalyzerInput+" --subset "+benchmarkSubset_Path
        postAnalysisResults_Path = getAllAfterLastChar(runBashCommand(bashCommand), "postAnalysis finished.").strip()
    else:
        if not args.excludeBasedonDistancePath:
            bashCommand = pythonPath+" RunPostAnalyzer.py --pythonPath "+pythonPath+" "+runPostAnalyzerInput+" --subset "+benchmarkSubset_Path
        else:
            bashCommand = pythonPath+" RunPostAnalyzer.py --pythonPath "+pythonPath+" "+runPostAnalyzerInput+" --subset "+benchmarkSubset_Path + " -excludeAfterDistance " + args.excludeAfterDistance + " --excludeBasedonDistancePath " + args.excludeBasedonDistancePath
        if excludefirstKModes:
            bashCommand = bashCommand + " --excludefirstKModes " + str(excludefirstKModes)
        if args.excludeFromList:
            bashCommand = bashCommand + " --excludeFromList " + args.excludeFromList
        postAnalysisResults_Path = getAllAfterLastChar(runBashCommand(bashCommand), "postAnalysis finished.").strip()
    print postAnalysisResults_Path

    # MakePlotfile.py
    if not args.dontCompare:
        bashCommand = pythonPath+" MakePlotFile.py --resultsDescriptor "+args.resultsDescriptor+" "+combinedIndividualResults_Path+" "+args.canonicalResults_Path
    else:
        # Do not compare results with optional (for instance canonical) results, just analyze and plot them
        bashCommand = pythonPath+" MakePlotFile.py --resultsDescriptor "+args.resultsDescriptor+" "+combinedIndividualResults_Path
    plotFile_Path = runBashCommand(bashCommand).strip()
    print plotFile_Path
    
    # AnalysisPlotter.py
    if args.plotOnlyReduction:
        bashCommand = pythonPath+" "+AnalysisPlotter_Path+" --plotOnlyReduction "+plotFile_Path+" "+postAnalysisResults_Path
    else:
        bashCommand = pythonPath+" "+AnalysisPlotter_Path+" "+plotFile_Path+" "+postAnalysisResults_Path
    plotFolderPath = getAllAfterLastChar(runBashCommand(bashCommand), "postPlot finished.").strip()
    if args.isComplex:    
        bashCommand = pythonPath+" "+AnalysisPlotter_Path+" --plotOnlyReduction --plotLRMS "+plotFile_Path+" "+postAnalysisResults_Path
        runBashCommand(bashCommand)
    print plotFolderPath
    
    if not args.dontCompare:
        # SuccessFailureAnalyzer.py        
        ### old bashCommand = pythonPath+" "+SuccessFailureAnalyzer_Path+" "+combinedIndividualResults_Path+ " "+args.canonicalResults_Path+" "+postAnalysisResults_Path
        if not args.excludeBasedonDistancePath:
            bashCommand = pythonPath+" "+SuccessFailureAnalyzer_Path+" -resultsPath "+combinedIndividualResults_Path+ " -resultsPath_Interface "+combinedIndividualResultsInterface_Path+" -canonicalResultsPath "+args.canonicalResults_Path+" -canonicalResultsPath_Interface "+args.canonicalResultsInterface_Path+" -outputPath "+postAnalysisResults_Path
        else:
            bashCommand = pythonPath+" "+SuccessFailureAnalyzer_Path+" -resultsPath "+combinedIndividualResults_Path+ " -resultsPath_Interface "+combinedIndividualResultsInterface_Path+" -canonicalResultsPath "+args.canonicalResults_Path+" -canonicalResultsPath_Interface "+args.canonicalResultsInterface_Path+" -outputPath "+postAnalysisResults_Path + " -excludeAfterDistance " + args.excludeAfterDistance + " --excludeBasedonDistancePath " + args.excludeBasedonDistancePath                                 
        successFailureAnalyzerOutput_Path = getAllAfterLastChar(runBashCommand(bashCommand), "successFailureAnalyzer finished.").strip()
        print successFailureAnalyzerOutput_Path
         
        # SuccessFailurePlotter.py
        if args.isComplex:
            bashCommand = pythonPath+" "+SuccessFailurePlotter_Path+" --isComplex --categoryMultiplier "+categoryMultiplier+" "+successFailureAnalyzerOutput_Path+" "+benchmarkSubset_Path+" "+taxonomyPath+" "+plotFolderPath
        else:
            bashCommand = pythonPath+" "+SuccessFailurePlotter_Path+" --categoryMultiplier "+categoryMultiplier+" "+successFailureAnalyzerOutput_Path+" "+sf_benchmarkSubsetPath+" "+taxonomyPath+" "+plotFolderPath
        runBashCommand(bashCommand)
        
#         #Plot AUC difference
#         bashCommand = pythonPath+" "+VarargsPlotter_Path+" --resultsDenoteAreas -title quality_measure_protein -xLabel modes -yLabel measure -outputFile "+postAnalysisResults_Path+"plots/whole_measure "+"--colorCode rgb "+postAnalysisResults_Path + "successiveAUCs/RMSD_auc_normalized_successiveProteinAll.txt " +postAnalysisResults_Path + "successiveAUCs/RMSD_auc_successiveProteinAll_stepPoints.txt " + "all " +postAnalysisResults_Path + "successiveAUCs/RMSD_auc_normalized_successiveProteinDifficult.txt " +postAnalysisResults_Path + "successiveAUCs/RMSD_auc_successiveProteinDifficult_stepPoints.txt " + "difficult "+postAnalysisResults_Path + "successiveAUCs/RMSD_auc_normalized_successiveProteinMedium.txt " +postAnalysisResults_Path + "successiveAUCs/RMSD_auc_successiveProteinMedium_stepPoints.txt " + "medium "+postAnalysisResults_Path + "successiveAUCs/RMSD_auc_normalized_successiveProteinRigid.txt " + postAnalysisResults_Path + "successiveAUCs/RMSD_auc_successiveProteinRigid_stepPoints.txt " + "rigid "
#         runBashCommand(bashCommand)
#         bashCommand = pythonPath+" "+VarargsPlotter_Path+" --resultsDenoteAreas -title quality_measure_interface -xLabel modes -yLabel measure -outputFile "+postAnalysisResults_Path+"plots/interface_measure "+"--colorCode rgb "+postAnalysisResults_Path + "successiveAUCs/RMSD_auc_normalized_successiveInterfaceAll.txt " +postAnalysisResults_Path + "successiveAUCs/RMSD_auc_successiveInterfaceAll_stepPoints.txt " + "all " +postAnalysisResults_Path + "successiveAUCs/RMSD_auc_normalized_successiveInterfaceDifficult.txt " +postAnalysisResults_Path + "successiveAUCs/RMSD_auc_successiveInterfaceDifficult_stepPoints.txt " + "difficult "+postAnalysisResults_Path + "successiveAUCs/RMSD_auc_normalized_successiveInterfaceMedium.txt " +postAnalysisResults_Path + "successiveAUCs/RMSD_auc_successiveInterfaceMedium_stepPoints.txt " + "medium "+postAnalysisResults_Path + "successiveAUCs/RMSD_auc_normalized_successiveInterfaceRigid.txt " + postAnalysisResults_Path + "successiveAUCs/RMSD_auc_successiveInterfaceRigid_stepPoints.txt " + "rigid "
#         runBashCommand(bashCommand)                
#         print "finished plotting AUC difference"
        
    # MakeLatexGallery.py
    if args.plotOnlyReduction:
        print "combining into latex gallery"
        bashCommand = pythonPath+" MakeLatexGallery.py "+plotFolderPath
        runBashCommand(bashCommand)

if __name__ == '__main__':
    main()