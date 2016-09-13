#!/bin/bash
#
# Visualize the output of Cprofile / pstats
# Required packages: gprof2dot, graphviz (for dot)
#
# Arguments:
#	$1: Cprofile / pstats output file
#
# Result:
#	$1.png file with a visualization of the Cprofile / pstats output file
gprof2dot -f pstats $1 | dot -Tpng -o $1.png