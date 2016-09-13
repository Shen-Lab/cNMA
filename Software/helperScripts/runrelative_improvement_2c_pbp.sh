#!/bin/bash

# paths
outputPrefix="/home/oliwa/data/output/individual/2c/alignment_beta/results_use_10AA/"
outputPath=$outputPrefix"results_use_relativeImprovement_final_05_pbp"

# 2x experiment, 2b or 2c
x2x="2c"
against="canonical"

#run quality measure plotter
resolution="bb"
model="2k"
customHIndividuals="submatrix HC_subvector HC_subvector_lambdaRTrue"
projectionString="projectionTrue"
align="beta"
projectionstyle="full"
gamma="0.5"

wholeOrInterfaceSet="whole interface"
difficulties="all difficult medium rigid"

for wholeOrInterface in $wholeOrInterfaceSet
do
  for difficulty in $difficulties
  do
    for customHIndividual in $customHIndividuals
    do
        if [ "$wholeOrInterface" == "whole" ]
        then
            wholeOrInterfaceAUCString="Whole"
            wholeOrInterfaceAUCStringProtein="Protein"
        else
            wholeOrInterfaceAUCString="Interface"
            wholeOrInterfaceAUCStringProtein="Interface"
        fi
      
        f6=$outputPrefix$x2x"_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_d6.0_gamma"$gamma"_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"_results/relativePerformanceAgainst_canonical/"
        f8=$outputPrefix$x2x"_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_d8.0_gamma"$gamma"_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"_results/relativePerformanceAgainst_canonical/"
        f10=$outputPrefix$x2x"_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_d10.0_gamma"$gamma"_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"_results/relativePerformanceAgainst_canonical/"
        f12=$outputPrefix$x2x"_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_d12.0_gamma"$gamma"_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"_results/relativePerformanceAgainst_canonical/"
        f15g05=$outputPrefix$x2x"_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_d15.0_gamma"$gamma"_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"_results/relativePerformanceAgainst_canonical/"
        f15g1=$outputPrefix$x2x"_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_d15.0_gamma1.0_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"_results/relativePerformanceAgainst_canonical/"

        y="RMSD_relativePerformance_"$wholeOrInterfaceAUCString""$difficulty"_pbp.txt"
        
        difficulty_uppercase="$(echo "$difficulty" | sed 's/.*/\u&/')"
        x="RMSD_relativePerformance_"$wholeOrInterfaceAUCStringProtein"Difficult_stepPoints.txt"
        
        d6="d_=_6.0_g="$gamma""
        d8="d_=_8.0_g="$gamma""
        d9="d_=_9.0_g="$gamma""
        d10="d_=_10.0_g="$gamma""
        d11="d_=_11.0_g="$gamma""
        d12="d_=_12.0_g="$gamma""
        d15g05="d_=_15.0_g="$gamma""
        d15g1="d_=_15.0_g=1.0"
        
        mkdir -p $outputPath"/"
        title=$x2x"_"$customHIndividual"_comparison_"$against
        # /home/oliwa/anaconda/bin/python /home/oliwa/workspace/TNMA1/src/VarargsPlotter.py -k 80 -ymin -100 -ymax 100 --colorCode lastTwoBlack -title $title"_"$wholeOrInterface"_"$difficulty -xLabel "Modes" -yLabel "Average Improvement (%)" -outputFile $outputPath"/"$title""$wholeOrInterface""$difficulty"all80_100" $f6$y $f6$x $d6 $f8$y $f8$x $d8 $f10$y $f10$x $d10 $f12$y $f12$x $d12 $f15g05$y $f15g05$x $d15g05 $f15g1$y $f15g1$x $d15g1
        /home/oliwa/anaconda/bin/python /home/oliwa/workspace/TNMA1/src/VarargsPlotterInset.py -k 80 -ymin -50 -ymax 50 --colorCode lastTwoBlack -title $title"_"$wholeOrInterface"_"$difficulty -xLabel "Modes" -yLabel "Average Improvement (%)" -outputFile $outputPath"/"$title""$wholeOrInterface""$difficulty"all80" $f6$y $f6$x $d6 $f8$y $f8$x $d8 $f10$y $f10$x $d10 $f12$y $f12$x $d12 $f15g05$y $f15g05$x $d15g05 $f15g1$y $f15g1$x $d15g1
    done
  done
done