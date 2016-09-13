#!/bin/bash

# paths
outputPrefix="/home/oliwa/data/output/complex/"

# parameters
# resolutionSet="bb all"
# models="2k"
# ds="6.0 8.0 10.0"
# gammas="0.5 1.0"
# projections="True False"

#run analysis and plotter on canonical 1k1k
# echo /home/oliwa/anaconda/bin/python Results2Plots.py --dontCompare --dontCombine --isComplex --pythonPath /home/oliwa/anaconda/bin/python --uniqueResultsPrefix_Path $outputPrefix"2c_complex_1k1k_whole_bb" --uniqueResultsPrefixInterface_Path $outputPrefix"2c_complex_1k1k_interface_bb" --resultsDescriptor Uk1k > runme.sge
# qsub -cwd runme.sge
# echo /home/oliwa/anaconda/bin/python Results2Plots.py --dontCompare --dontCombine --isComplex --pythonPath /home/oliwa/anaconda/bin/python --uniqueResultsPrefix_Path $outputPrefix"2c_complex_1k1k_whole_all" --uniqueResultsPrefixInterface_Path $outputPrefix"2c_complex_1k1k_interface_all" --resultsDescriptor Uk1k > runme.sge
# qsub -cwd runme.sge

# for model in $models
# do
#   for resolution in $resolutionSet
#   do
#     foldernameWhole=$outputPrefix"2c_complex_"$model"_whole_"$resolution"/"
#     foldernameInterface=$outputPrefix"2c_complex_"$model"_interface_"$resolution"/"
#     compareAgainstfoldernameWhole=$outputPrefix"2c_complex_1k1k_whole_"$resolution"/"
#     compareAgainstfoldernameInterface=$outputPrefix"2c_complex_1k1k_interface_"$resolution"/"
#     
#     #run analysis and plotter
#     echo /home/oliwa/anaconda/bin/python Results2Plots.py --dontCombine --isComplex --pythonPath /home/oliwa/anaconda/bin/python --uniqueResultsPrefix_Path $foldernameWhole --uniqueResultsPrefixInterface_Path $foldernameInterface --canonicalResults_Path $compareAgainstfoldernameWhole --canonicalResultsInterface_Path $compareAgainstfoldernameInterface --resultsDescriptor U$model > runme.sge
#     qsub -cwd runme.sge
#     " " > runme.sge 
#   done
# done

resolutionSet="bb"
models="2k"
ds="6.0"
gammas="0.5"
projections="True False"

for model in $models
do
  for resolution in $resolutionSet
  do
    for d in $ds
    do
      for gamma in $gammas
      do
        for projection in $projections
        do
          if [ "$projection" == "True" ]
          then
            projectionString="_projectionTrue"
          else
            projectionString=""
          fi
          
          foldernameWhole=$outputPrefix"2c_complex_"$model"_whole_"$resolution"_HC_U1_d"$d"_gamma"$gamma""$projectionString"/"
          foldernameInterface=$outputPrefix"2c_complex_"$model"_interface_"$resolution"_HC_U1_d"$d"_gamma"$gamma""$projectionString"/"
          compareAgainstfoldernameWhole=$outputPrefix"2c_complex_1k1k_whole_"$resolution"/"
          compareAgainstfoldernameInterface=$outputPrefix"2c_complex_1k1k_interface_"$resolution"/"
          
          #run analysis and plotter
          echo /home/oliwa/anaconda/bin/python Results2Plots.py --dontCombine --isComplex --pythonPath /home/oliwa/anaconda/bin/python --uniqueResultsPrefix_Path $foldernameWhole --uniqueResultsPrefixInterface_Path $foldernameInterface --canonicalResults_Path $compareAgainstfoldernameWhole --canonicalResultsInterface_Path $compareAgainstfoldernameInterface --resultsDescriptor U$model"d"$d"k"$gamma"p"$projection > runme.sge
          qsub -cwd runme.sge
        done
      done
    done
  done
done