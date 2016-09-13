#!/bin/bash

# paths
outputPrefix="/home/oliwa/data/output/individual/2b/"

# parameters
alignstyle="align2bindividual"

# run analysis and plotter on canonical 1k1k
/home/oliwa/anaconda/bin/python Results2Plots.py --dontCompare --dontCombine --plotOnlyReduction --pythonPath /home/oliwa/anaconda/bin/python --uniqueResultsPrefix_Path $outputPrefix"2b_individual_canonical_whole_bb_"$alignstyle --uniqueResultsPrefixInterface_Path $outputPrefix"2b_individual_canonical_interface_bb_"$alignstyle --resultsDescriptor Uk1k

