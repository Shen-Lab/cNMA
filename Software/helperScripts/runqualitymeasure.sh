#!/bin/bash

# paths
outputPrefix="/home/oliwa/data/output/individual/2c/alignment_beta/"

# 2x experiment, 2b or 2c
x2x="2c"

#run quality measure plotter
resolution="bb"
model="2k"
customHIndividual="submatrix"
projectionString="projectionTrue"
align="beta"
projectionstyle="full"

title=$x2x"_"$customHIndividual"_comparison"
wholeOrInterfaceSet="whole interface"
difficulties="All Difficult Medium Rigid"

for wholeOrInterface in $wholeOrInterfaceSet
do
  for difficulty in $difficulties
  do
    if [ "$wholeOrInterface" == "whole" ]
    then
        wholeOrInterfaceAUCString="Protein"
    else
        wholeOrInterfaceAUCString="Interface"
    fi
  
    f6=$outputPrefix$x2x"_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_d6.0_gamma0.5_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"_results/successiveAUCs/"
    f8=$outputPrefix$x2x"_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_d8.0_gamma0.5_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"_results/successiveAUCs/"
    f9=$outputPrefix$x2x"_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_d9.0_gamma0.5_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"_results/successiveAUCs/"
    f10=$outputPrefix$x2x"_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_d10.0_gamma0.5_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"_results/successiveAUCs/"
    #f11=$outputPrefix$x2x"_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_d11.0_gamma0.5_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"_results/successiveAUCs/"
    f12=$outputPrefix$x2x"_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_d12.0_gamma0.5_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"_results/successiveAUCs/"
    f15=$outputPrefix$x2x"_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_d15.0_gamma1.0_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"_results/successiveAUCs/"
    
    y="RMSD_auc_normalized_successive"$wholeOrInterfaceAUCString""$difficulty".txt"
    x="RMSD_auc_successive"$wholeOrInterfaceAUCString""$difficulty"_stepPoints.txt"
    
    d6="d_=_6.0"
    d8="d_=_8.0"
    d9="d_=_9.0"
    d10="d_=_10.0"
    #d11="d_=_11.0"
    d12="d_=_12.0"
    d15="d_=_15.0"
    
    /home/oliwa/anaconda/bin/python /home/oliwa/workspace/TNMA1/src/VarargsPlotter.py -title $title"_"$wholeOrInterface"_"$difficulty -xLabel "modes" -yLabel "measure" -outputFile $outputPrefix"results_qualitymeasures/"$title""$wholeOrInterface""$difficulty $f6$y $f6$x $d6 $f8$y $f8$x $d8 $f9$y $f9$x $d9 $f10$y $f10$x $d10 $f12$y $f12$x $d12 $f15$y $f15$x $d15
  done
done