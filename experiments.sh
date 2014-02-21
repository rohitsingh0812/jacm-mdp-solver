#!/bin/bash
if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters
Usage:
bash experiemnts.sh <low> <high>
e.g.
bash experiemnts.sh 3 10
"
fi
low=$1
up=$2

pall=0.25 #probability of a process being hungry (same for all)

#the sequential probabilities start from 0.05 and add 0.07 for each next one
parr=( 0.05 0.12 0.19 0.26 0.33 0.40 0.47 0.54 0.61 0.68 )
for i in `seq $low $up`;
do
    bash run_gcc.sh $i 1 $pall
    bash run_gcc.sh $i 2 "${parr[@]:0:$i}"
done    
