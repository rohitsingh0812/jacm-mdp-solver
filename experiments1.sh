#!/bin/bash
if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters
Usage:
bash experiments1.sh <low> <high> 
e.g.
bash experiments1.sh 2 8
"
fi
low=$1
up=$2
echo "" > tables1_true.txt
echo "" > tables1_false.txt

for i in `seq $low $up`;
do
	echo "run_exp1.sh $i true"
    bash run_exp1.sh $i "true"
    echo "run_exp1.sh $i false"    
    bash run_exp1.sh $i "false"
done    
