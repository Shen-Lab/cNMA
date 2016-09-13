#!/bin/bash
# Run 2c for the complex
# Arguments:
#   $1: the input file with: prot1a prot2a prot1b prot2b

# Paths
pythonPath="/home/oliwa/anaconda/bin/python"
NMAPath="/home/oliwa/workspace/TNMA1/src/"
NMAProgram="NMAUnified.py --profile"
inputFile=$1
configPath="/home/oliwa/workspace/TNMA1/src/TNMAConfigurations/runscripts_Individual_NMAUnified/"
configFile="Configurations_NMAUnifiedIndividualTemplate.py"
tempPath="/tmp/oliwa_\$JOB_ID"

# Create batch indexes of the input
inputLines=$(wc $inputFile | awk {'print $1'})
batchSize=$(($inputLines/34))

# parameters
resolutionSet="bb"
measuresOnWhole="True False"
complexRMSDreduction="2k"
potential="HC_U1"
whichCustomHIndividual="HC_subvector"
distanceCutoffs="15.0"
gammas="0.25 0.5 1.0"
projections="True"
alignstyle="beta"
projectionstyles="full"
rescaleEigenvalues="True"

for resolution in $resolutionSet
do
  for measureOnWhole in $measuresOnWhole
  do
    for RMSDreduction in $complexRMSDreduction
    do
      for d in $distanceCutoffs
      do
        for gamma in $gammas
        do
          for projection in $projections
          do
            for projectionstyle in $projectionstyles
            do
              for customHIndividual in $whichCustomHIndividual
              do
                for lambdaR in $rescaleEigenvalues
                do
                    ###customConfiguration="2c_complex_bb_whole_1k1k"
                    if [ $measureOnWhole = "True" ];
                    then
                    measureOnWholeLiteral="whole"
                    else
                    measureOnWholeLiteral="interface"
                    fi      
                    # create the new configuration name
                    experimentNamePrefix="2c_individual_"$RMSDreduction"_"$measureOnWholeLiteral"_"$resolution"_"$potential"_model"$customHIndividual"_lambdaR"$lambdaR"_d"$d"_gamma"$gamma"_projection"$projection"_projectionStyle"$projectionstyle"_align"$alignstyle
                    
                    # fill the configuration template
                    cat $configPath$configFile | 
                    sed "s/self.experimentNamePrefix = \"x\"/self.experimentNamePrefix = \"$experimentNamePrefix\"/g" | 
                    sed "s/self.measuresOnWhole = x/self.measuresOnWhole = $measureOnWhole/g" |
                    sed "s/self.align = x/self.align = \"$alignstyle\"/g" |
                    sed "s/self.whatAtomsToMatch = \"x\"/self.whatAtomsToMatch = \"$resolution\"/g" |
                    sed "s/self.complexRMSDreduction = \"x\"/self.complexRMSDreduction = \"$customHIndividual\"/g" |
                    sed "s/self.whichCustomHIndividual = \"x\"/self.whichCustomHIndividual = \"$customHIndividual\"/g" |
                    sed "s/self.whichCustomHC = \"x\"/self.whichCustomHC = \"HC_U1\"/g" |
                    sed "s/self.customHRdistance = 10.0/self.customHRdistance = $d/g" |
                    sed "s/self.customForceConstant = 0.5/self.customForceConstant = $gamma/g" |
                    sed "s/self.projectHessian = False/self.projectHessian = $projection/g" |
                    sed "s/self.rescaleEigenvalues = False/self.rescaleEigenvalues = $lambdaR/g" | 
                    sed "s/self.projectionStyle = \"x\"/self.projectionStyle = \"$projectionstyle\"/g" > $configPath$experimentNamePrefix.p2

                    if [ $RMSDreduction = "2k" ];
                    then
                    cat $configPath$experimentNamePrefix.p2 |
                    sed "s/self.customH = x/self.customH = True/g" > $configPath$experimentNamePrefix.py
                    else
                    cat $configPath$experimentNamePrefix.p2 |
                    sed "s/self.customH = x/self.customH = False/g" > $configPath$experimentNamePrefix.py
                    fi 
                    rm $configPath$experimentNamePrefix.p2
                    
                    echo $configPath$experimentNamePrefix
                    rm runme.sge
                    rm proteinBatch.txt
                    rm proteinBatchNoPath.txt
                    
                    lineCounter=0
                    while read line; do
                        lineNoPaths=$(echo $line | sed "s,/home/oliwa/nmaData/2c/,,g" | sed "s,/home/oliwa/nmaData/benchmark/,,g")
                        lineCounter=$((lineCounter+1))
                        echo $line >> proteinBatch.txt
                        echo $lineNoPaths >> proteinBatchNoPath.txt
                        currentBatchItem=$((lineCounter % batchSize))
                        echo $lineNoPaths | sed "s,^,$pythonPath $NMAProgram --outputPath "." $experimentNamePrefix.py ,g" >> runme.sge
                        if (($currentBatchItem==0)); # can be done to run 2 lines of the input with (($lineCounter==2)), else use (($currentBatchItem==0)) for batchSize
                        then
                        echo $lineCounter
                        # qsub customize the batchscript
                        echo "#$ -S /bin/bash" > runme.sge
                        echo "#$ -e /tmp/" >> runme.sge
                        echo "#$ -o /tmp/" >> runme.sge                   
                        echo rm -rf $tempPath >> runme.sge
                        echo mkdir -p $tempPath >> runme.sge
                        echo cp $NMAPath*.py $tempPath >> runme.sge
                        echo cp -r $NMAPath\helperScripts $tempPath >> runme.sge
                        echo cp $configPath$experimentNamePrefix.py $tempPath >> runme.sge
                        awk '{ for(count=1; count<=NF; count++) print $count }' proteinBatch.txt > proteinPaths.txt
                        for element in `<proteinPaths.txt`; do echo cp $element $tempPath >> runme.sge ; done;
                        echo cd $tempPath >> runme.sge
                        cat proteinBatchNoPath.txt | sed "s,^,$pythonPath $NMAProgram --outputPath "." $experimentNamePrefix.py ,g" >> runme.sge
                        echo echo JOB_ID: \$JOB_ID >> runme.sge
                        echo echo JOB_NAME: \$JOB_NAME >> runme.sge
                        echo echo batchline: $lineCounter
                        # copy results to real output folder
                        realOutputPath=$(cat $configPath$experimentNamePrefix.py | grep "self.outputPath" | sed "s/        self.outputPath = //g" | sed "s/\"//g")
                        echo cp -r $experimentNamePrefix $realOutputPath >> runme.sge
                        # clean the o and e files and copy them to the real output folder
                        echo cd /tmp/ >> runme.sge
                        echo 'errorfileName=$(ls *.e$JOB_ID)' >> runme.sge
                        echo 'cat $errorfileName > temp_e$JOB_ID\.txt ; cat temp_e$JOB_ID\.txt | grep -v "@> WARNING failed to parse beta-factor" > fixed_e$JOB_ID\.txt' >> runme.sge
                        echo 'mv fixed_e$JOB_ID\.txt $errorfileName' >> runme.sge
                        echo 'rm temp_e$JOB_ID\.txt' >> runme.sge
                        echo cp *.e\$JOB_ID $realOutputPath\sgefiles/ >> runme.sge
                        echo cp *.o\$JOB_ID $realOutputPath\sgefiles/ >> runme.sge
                        # cleanup 
                        echo cd /tmp/ >> runme.sge
                        echo rm -rf $tempPath >> runme.sge
                        echo rm *.e\$JOB_ID >> runme.sge
                        echo rm *.o\$JOB_ID >> runme.sge
                        # run the configuration
                        cp runme.sge run$experimentNamePrefix.sge
                        qsub -cwd -N "j2c_individual_"$customHIndividual"_"$measureOnWholeLiteral"_"$resolution"_"$potential"_lambdaR"$lambdaR"_d"$d"_gamma"$gamma"_projection"$projection"_projectionStyle"$projectionstyle"_align"$alignstyle runme.sge
                        rm runme.sge
                        rm proteinPaths.txt
                        rm proteinBatch.txt
                        rm proteinBatchNoPath.txt
                        fi
                    done < $1
                done
              done
            done
          done
        done
      done
    done
  done
done