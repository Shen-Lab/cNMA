#!/bin/bash

# paths
outputPrefix="/home/oliwa/data/output/complex/"

# parameters
resolutionSet="bb"
models="2k"
alignstyle="align2bindividual"

# run analysis and plotter on canonical 1k1k
#/home/oliwa/anaconda/bin/python Results2Plots.py --dontCompare --dontCombine --isComplex --pythonPath /home/oliwa/anaconda/bin/python --uniqueResultsPrefix_Path $outputPrefix"2b_complex_1k1k_whole_bb_HC_0_"$alignstyle --uniqueResultsPrefixInterface_Path $outputPrefix"2b_complex_1k1k_interface_bb_HC_0_"$alignstyle --resultsDescriptor Uk1k
#/home/oliwa/anaconda/bin/python Results2Plots.py --dontCompare --dontCombine --isComplex --pythonPath /home/oliwa/anaconda/bin/python --uniqueResultsPrefix_Path $outputPrefix"2b_complex_1k1k_whole_all_HC_0_"$alignstyle --uniqueResultsPrefixInterface_Path $outputPrefix"2b_complex_1k1k_interface_all_HC_0_"$alignstyle --resultsDescriptor Uk1k

for model in $models
do
  for resolution in $resolutionSet
  do
    for align in $alignstyle
    do
        foldernameWhole=$outputPrefix"2b_complex_"$model"_whole_"$resolution"_HC_0_"$align"/"
        foldernameInterface=$outputPrefix"2b_complex_"$model"_interface_"$resolution"_HC_0_"$align"/"
        compareAgainstfoldernameWhole=$outputPrefix"2b_complex_1k1k_whole_"$resolution"_HC_0_"$align"/"
        compareAgainstfoldernameInterface=$outputPrefix"2b_complex_1k1k_interface_"$resolution"_HC_0_"$align"/"
        
        #run analysis and plotter
        /home/oliwa/anaconda/bin/python Results2Plots.py --dontCombine --isComplex --pythonPath /home/oliwa/anaconda/bin/python --uniqueResultsPrefix_Path $foldernameWhole --uniqueResultsPrefixInterface_Path $foldernameInterface --canonicalResults_Path $compareAgainstfoldernameWhole --canonicalResultsInterface_Path $compareAgainstfoldernameInterface --resultsDescriptor U$model
    done
  done
done