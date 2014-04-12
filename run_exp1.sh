#!/bin/bash
if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters
Usage:
bash run_exp1.sh <num_elevators> <weighted-rewards:true|false>
e.g.
bash run_exp1.sh 2 true
";
    exit 1;
fi
NF=$1

echo "#define NF $1
#define WEIGHTED_REWARDS $2
" > inputs1.h


DN="e1_$NF_$2" 
mkdir -p $DN

g++ -o p1 MDPex1.cpp

./p1 > ./$DN/output.txt

cat ./$DN/output.txt | tail -n 1 >> tables1_$2.txt


