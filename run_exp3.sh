#!/bin/bash
if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters
Usage:
bash run_exp3.sh <num_products>
e.g.
bash run_exp3.sh 2 
";
    exit 1;
fi
NF=$1
STRUP=`seq -s, 1 $NF`
STRDOWN=`seq -s, $NF -1 1`

echo "#define NF $1
double pall = 0.1; //MDPTYPE 1
int rewards[NF] = {$STRDOWN};
int steps[NF] = {$STRUP};
" > inputs3_temp.h

echo "class sysvars{
public:
  int prod;
  int step;
  bool repair;  
  int reward;
  sysvars(int p,int s, bool r, int rew){
    prod = p; step= s; repair = r; reward = rew;
  }
};

#define MDPTYPE 2
double pdiff[2] = {0.05,0.12}; //MDPTYPE 2
//#define NK 3
" >> inputs3_temp.h

mv inputs3_temp.h inputs3.h

DN="e3_$NF" 
mkdir -p $DN

g++ -o p3 MDPex3.cpp

./p3 > ./$DN/output.txt

cat ./$DN/output.txt | tail -n 1 >> tables3.txt


