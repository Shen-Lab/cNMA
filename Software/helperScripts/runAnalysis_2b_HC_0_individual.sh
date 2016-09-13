#!/bin/bash

# paths
outputPrefix="/home/oliwa/data/output/individual/2b/"

# parameters
alignstyle="align2bindividual"

# run analysis and plotter on canonical 1k1k
#/home/oliwa/anaconda/bin/python Results2Plots.py --dontCompare --dontCombine --pythonPath /home/oliwa/anaconda/bin/python --uniqueResultsPrefix_Path $outputPrefix"2b_individual_canonical_whole_bb_"$alignstyle --uniqueResultsPrefixInterface_Path $outputPrefix"2b_individual_canonical_interface_bb_"$alignstyle --resultsDescriptor Uk1k

#run analysis and plotter on U1 
resolutionSet="bb"
models="2k"
ds="17.0"
gammas="0.5"
whichCustomHIndividual="HC_subvector"
projections="True"
alignstyle="2bindividual"
projectionstyles="full"

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
          for align in $alignstyle
          do
            for customHIndividual in $whichCustomHIndividual
              do
                for projectionstyle in $projectionstyles
                do
                    if [ "$projection" == "True" ]
                    then
                        projectionString="projectionTrue"
                    else
                        projectionString="projectionFalse"
                    fi
                    #_projectionStylefull_alignmentL
                    
                    #2b_individual_2k_whole_bb_HC_U1_modelHC_subvector_lambdaRTrue_d12.0_gamma0.5_projectionTrue_projectionStylefull_align2bindividual
                    foldernameWhole=$outputPrefix"2b_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_lambdaRTrue_d"$d"_gamma"$gamma"_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"/"
                    foldernameInterface=$outputPrefix"2b_individual_"$model"_interface_"$resolution"_HC_U1_model"$customHIndividual"_lambdaRTrue_d"$d"_gamma"$gamma"_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"/"
                    compareAgainstfoldernameWhole=$outputPrefix"2b_individual_canonical_whole_bb_align2bindividual/"
                    compareAgainstfoldernameInterface=$outputPrefix"2b_individual_canonical_interface_bb_align2bindividual/"
                    
                    #run analysis and plotter
                    /home/oliwa/anaconda/bin/python Results2Plots.py --dontCombine --categoryMultiplier 1.0 --pythonPath /home/oliwa/anaconda/bin/python --uniqueResultsPrefix_Path $foldernameWhole --uniqueResultsPrefixInterface_Path $foldernameInterface --canonicalResults_Path $compareAgainstfoldernameWhole --canonicalResultsInterface_Path $compareAgainstfoldernameInterface --resultsDescriptor $customHIndividual"d"$d"k"$gamma"p"$projection
                    
              done
            done
          done
        done
      done
    done
  done
done


                    #########compareAgainstfoldernameWhole=$outputPrefix"2b_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_lambdaRTrue_d15.0_gamma1.0_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"/"
                    #########compareAgainstfoldernameInterface=$outputPrefix"2b_individual_"$model"_interface_"$resolution"_HC_U1_model"$customHIndividual"_lambdaRTrue_d15.0_gamma1.0_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"/"  