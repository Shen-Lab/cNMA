#!/bin/bash 
# run the Results2Plots.py

# Paths
pythonPath="/home/oliwa/anaconda/bin/python"
NMAPath="/home/oliwa/workspace/TNMA1/src/helperScripts/Results2Plots.py"

# --dontCombine --dontCompare --isComplex --uniqueResultsPrefix_Path /home/oliwa/nmaData/output/testing/2b_all_Complex_kx_Ax_lambdaC_HC_1k1k/ --uniqueResultsPrefixInterface_Path /home/oliwa/nmaData/output/testing/2b_all_Complex_kx_Ax_lambdaC_HC_1k1k_interface/ --resultsDescriptor U1test
# --dontCombine --dontCompare --isComplex --uniqueResultsPrefix_Path /home/oliwa/data/test/2c_complex_1k1k_whole_bb/ --uniqueResultsPrefixInterface_Path /home/oliwa/data/test/2c_complex_1k1k_interface_bb/ --resultsDescriptor U1test

# run the configuration
# canonical
# $pythonPath $NMAPath --dontCombine --dontCompare --pythonPath /home/oliwa/anaconda/bin/python --uniqueResultsPrefix_Path /home/oliwa/nmaData/data/output/individual/2b_individual_canonical_whole_bb_align2bindividual/ --uniqueResultsPrefixInterface_Path /home/oliwa/nmaData/data/output/individual/2b_individual_canonical_interface_bb_align2bindividual/ --resultsDescriptor U1test > runme.sge
$pythonPath $NMAPath --dontCombine --pythonPath /home/oliwa/anaconda/bin/python --uniqueResultsPrefix_Path /home/oliwa/nmaData/data/output/individual/2b_individual_2k_whole_bb_HC_U1_modelHC_subvector_d8.0_gamma0.5_projectionFalse_projectionStylefull_align2bindividual --uniqueResultsPrefixInterface_Path /home/oliwa/nmaData/data/output/individual/2b_individual_2k_interface_bb_HC_U1_modelsubmatrix_d9.0_gamma0.5_projectionTrue_projectionStylefull_align2bindividual --canonicalResults_Path /home/oliwa/nmaData/data/output/individual/2b_individual_canonical_whole_bb_align2bindividual/ --canonicalResultsInterface_Path /home/oliwa/nmaData/data/output/individual/2b_individual_canonical_interface_bb_align2bindividual/ --resultsDescriptor HRtilde_d9_g05_pFull
# qsub -cwd runme.sge
# echo "" > runme.sge
    