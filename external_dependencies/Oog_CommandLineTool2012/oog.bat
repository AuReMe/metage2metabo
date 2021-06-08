REM Simple Oog batch script for windows/dos
REM (c) Matthias Reimann, 2012;
REM

java -Xmx1024M -jar Oog.jar -convert -inputfiles=edg -output=bbl -v -writeLog=0 -sortbyfilesize -checkequality -minsim=0.0 -min_max=1 -useGlobalClustering=1 -useGreedyClustering=1 -lat=0

REM Use the following line for help
REM java -Xmx1024M -jar Oog.jar -?
pause

