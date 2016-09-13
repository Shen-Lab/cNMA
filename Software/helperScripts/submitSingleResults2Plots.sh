#!/bin/bash 
# run the Results2Plots.py

# Paths
pythonPath="/home/oliwa/anaconda/bin/python"
NMAPath="/home/oliwa/workspace/TNMA1/src/helperScripts/Results2Plots.py"

# --dontCombine --dontCompare --isComplex --uniqueResultsPrefix_Path /home/oliwa/nmaData/output/testing/2b_all_Complex_kx_Ax_lambdaC_HC_1k1k/ --uniqueResultsPrefixInterface_Path /home/oliwa/nmaData/output/testing/2b_all_Complex_kx_Ax_lambdaC_HC_1k1k_interface/ --resultsDescriptor U1test
# --dontCombine --dontCompare --isComplex --uniqueResultsPrefix_Path /home/oliwa/data/test/2c_complex_1k1k_whole_bb/ --uniqueResultsPrefixInterface_Path /home/oliwa/data/test/2c_complex_1k1k_interface_bb/ --resultsDescriptor U1test

# run the configuration
# canonical
# echo -e "#!/bin/bash\n"$pythonPath $NMAPath --dontCombine --dontCompare --isComplex --pythonPath /home/oliwa/anaconda/bin/python --uniqueResultsPrefix_Path /home/oliwa/data/output/complex/2b_complex_1k1k_whole_bb_HC_0_align2bcomplex/ --uniqueResultsPrefixInterface_Path /home/oliwa/data/output/complex/2b_complex_1k1k_interface_bb_HC_0_align2bcomplex/ --resultsDescriptor U1test > runme.sge
echo -e "#!/bin/bash\n"$pythonPath $NMAPath --dontCombine --isComplex --pythonPath /home/oliwa/anaconda/bin/python --uniqueResultsPrefix_Path /home/oliwa/data/output/complex/2b_complex_2k_whole_bb_HC_0_align2bcomplex/ --uniqueResultsPrefixInterface_Path /home/oliwa/data/output/complex/2b_complex_2k_interface_bb_HC_0_align2bcomplex/ --canonicalResults_Path /home/oliwa/data/output/complex/2b_complex_1k1k_whole_bb_HC_0_align2bcomplex/ --canonicalResultsInterface_Path /home/oliwa/data/output/complex/2b_complex_1k1k_interface_bb_HC_0_align2bcomplex/ --resultsDescriptor HC0 > runme.sge
qsub -cwd runme.sge
# echo "" > runme.sge
    