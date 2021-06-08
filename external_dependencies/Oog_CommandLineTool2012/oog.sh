#!/bin/sh
#
# Run Oog from a jar file
# this is a linux-only version
#-------------------------------------------------------------------------------
java -Xmx12000M -jar Oog.jar -convert -inputfiles=edg -output=bbl -v -writeLog=0 -sortbyfilesize -checkequality -minsim=0.0 -min_max=1 -useGlobalClustering=1 -useGreedyClustering=1 -lat=0

# Use the following line for help
#java -Xmx1024M -jar Oog.jar -?
