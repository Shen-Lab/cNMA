#!/bin/bash

# paths
outputPrefix="/home/oliwa/data/output/individual/2c/alignment_beta/"

# parameters
alignstyle="alignbeta"

# run analysis and plotter on canonical 1k1k
#/home/oliwa/anaconda/bin/python Results2Plots.py --dontCompare --dontCombine --pythonPath /home/oliwa/anaconda/bin/python --uniqueResultsPrefix_Path $outputPrefix"2c_individual_canonical_whole_bb_HC_U1_modelcanonical_dx_gammax_projectionFalse_projectionStylefull_"$alignstyle --uniqueResultsPrefixInterface_Path $outputPrefix"2c_individual_canonical_interface_bb_HC_U1_modelcanonical_dx_gammax_projectionFalse_projectionStylefull_"$alignstyle --resultsDescriptor Ucanon

#run analysis and plotter on U1 
resolutionSet="bb"
models="2k"
ds="6.0 8.0 9.0 10.0 12.0"
gammas="0.5"
whichCustomHIndividual="HC_subvector"
projections="True"
alignstyle="beta"
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
                    foldernameWhole=$outputPrefix"2c_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_lambdaRTrue_d"$d"_gamma"$gamma"_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"/"
                    foldernameInterface=$outputPrefix"2c_individual_"$model"_interface_"$resolution"_HC_U1_model"$customHIndividual"_lambdaRTrue_d"$d"_gamma"$gamma"_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"/"
                    ##compareAgainstfoldernameWhole=$outputPrefix"2c_individual_canonical_whole_"$resolution"_HC_U1_modelcanonical_dx_gammax_projectionFalse_projectionStylefull_align"$align"/"
                    ##compareAgainstfoldernameInterface=$outputPrefix"2c_individual_canonical_interface_"$resolution"_HC_U1_modelcanonical_dx_gammax_projectionFalse_projectionStylefull_align"$align"/"
                    compareAgainstfoldernameWhole=$outputPrefix"2c_individual_"$model"_whole_"$resolution"_HC_U1_model"$customHIndividual"_lambdaRTrue_d15.0_gamma1.0_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"/"
                    compareAgainstfoldernameInterface=$outputPrefix"2c_individual_"$model"_interface_"$resolution"_HC_U1_model"$customHIndividual"_lambdaRTrue_d15.0_gamma1.0_"$projectionString"_projectionStyle"$projectionstyle"_align"$align"/"         
                    
                    #run analysis and plotter
                    echo /home/oliwa/anaconda/bin/python Results2Plots.py --dontCombine --categoryMultiplier 10.0 --pythonPath /home/oliwa/anaconda/bin/python --uniqueResultsPrefix_Path $foldernameWhole --uniqueResultsPrefixInterface_Path $foldernameInterface --canonicalResults_Path $compareAgainstfoldernameWhole --canonicalResultsInterface_Path $compareAgainstfoldernameInterface --resultsDescriptor $customHIndividual"d"$d"k"$gamma"p"$projection"ad" > runme.sge
                    qsub -cwd runme.sge
              done
            done
          done
        done
      done
    done
  done
done