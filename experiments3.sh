#!/bin/bash
if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters
Usage:
bash experiments3.sh <low> <high> 
e.g.
bash experiments3.sh 2 12
"
fi
low=$1
up=$2
echo "" > tables3.txt

for i in `seq $low $up`;
do
	echo "run_exp3.sh $i"
    bash run_exp3.sh $i
done    
