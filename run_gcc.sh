#!/bin/bash

NF=$1 #3-15
MDPTYPE=$2 #1 or 2
	echo "#define NF $NF
#define MDPTYPE $MDPTYPE" > ./temp.h

pall="0.05"

if [ "$MDPTYPE" = "1" ];
then
	pall=$3 #double value like 0.05
	echo "double pall = $pall; //MDPTYPE 1
double pdiff[NF];" >> ./temp.h
else
	echo -n "double pall = $pall; //MDPTYPE 1
double pdiff[NF] = {" >> ./temp.h

	for i in ${@:3} #${@:2} is the same array less the first element.
	do
		echo -n $i"," >> ./temp.h
	done
	echo -n "}; //MDPTYPE 2
" >> ./temp.h

	#followed by $NF values for probabilities if MDPTYPE is 2	
fi

cat temp.h | sed 's/,}/}/g' > inputs.h


g++ -o runner MDPex2.cpp

fldr=${NF}_${MDPTYPE}

mkdir -p $fldr

./runner > ./$fldr/output.txt

cd $fldr
for i in `ls strategy*.dot | sed 's/.dot//g'`; 
do
	dot -Tpng $i.dot > $i.png
done
