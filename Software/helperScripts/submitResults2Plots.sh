#!/bin/bash 
# run the Results2Plots.py for different configuration files, which are created based on a template

# Paths
pythonPath="/home/oliwa/anaconda/bin/python"
NMAPath="/home/oliwa/workspace/TNMA1/src/helperScripts/Results2Plots.py"
examplePDB="/home/oliwa/nmaData/benchmark/"
configPath="/home/oliwa/workspace/TNMA1/src/TNMAConfigurations/2B_Complex_1k1k_HC0_HC_U1/"
configFile_Template="2b_X_Complex_kX_AX_lambdaX_HCX.py"
outputPath="/home/oliwa/nmaData/output/300modes/"

# parameters
resolutionSet="bb all"
lambdaSet="C"
kSet="025 05"
dSet="6 8 10 12 15"

for resolution in $resolutionSet
do
  for lambda in $lambdaSet
  do
    # run the 1k1k configuration
    # create the new configuration name
    customConfigurationCanonical=$(echo $configFile_Template | \
    sed "s/2b_X/2b_$resolution/1" | \
    sed "s/_kX/_kx/1" | \
    sed "s/_AX/_Ax/1" | \
    sed "s/_lambdaX/_lambda$lambda/1" | \
    sed "s/_HCX/_HC_1k1k/1")
    customConfigurationCanonicalNoExtension=$(echo $customConfigurationCanonical | sed 's/...$//')
    
    # run the configuration
    # echo -e "#!/bin/bash\n"$pythonPath $NMAPath --dontCompare --isComplex --pythonPath $outputPath$customConfigurationCanonicalNoExtension 1k1k > runme.sge
    # qsub -cwd runme.sge
    # echo "" > runme.sge
    
    # run the HC_0 configuration (order effect)
    # create the new configuration name
    customConfiguration=$(echo $configFile_Template | \
    sed "s/2b_X/2b_$resolution/1" | \
    sed "s/_kX/_kx/1" | \
    sed "s/_AX/_Ax/1" | \
    sed "s/_lambdaX/_lambda$lambda/1" | \
    sed "s/_HCX/_HC_0/1")
    customConfigurationNoExtension=$(echo $customConfiguration | sed 's/...$//')
    
    # run the configuration
    # echo $pythonPath $NMAPath --isComplex --pythonPath $outputPath$customConfigurationNoExtension HC0 $outputPath$customConfigurationCanonicalNoExtension > runme.sge
    # qsub -cwd runme.sge
    # echo "" > runme.sge
    
    # run HC_U1 configurations
    for k in $kSet
    do
      for d in $dSet
	do
	  # create the new configuration name
	  customConfiguration=$(echo $configFile_Template | \
	  sed "s/2b_X/2b_$resolution/1" | \
	  sed "s/_kX/_k$k/1" | \
	  sed "s/_AX/_A$d/1" | \
	  sed "s/_lambdaX/_lambda$lambda/1" | \
	  sed "s/_HCX/_HC_U1/1")
	  customConfigurationNoExtension=$(echo $customConfiguration | sed 's/...$//')
	  
	  # run the configuration
	  echo $pythonPath $NMAPath --isComplex --pythonPath $outputPath$customConfigurationNoExtension U1_kd$k$d $outputPath$customConfigurationCanonicalNoExtension > runme.sge
	  qsub -cwd runme.sge
	  echo "" > runme.sge
      done
    done
  done
done
