#!/bin/bash


if [ $# -ne 1 ]; then
echo ERROR in section.sh: one argument is needed, e.g.: section.sh 5 cross sects the file 00005.hc
exit
fi

var="density"
px=0
py=0
qx=0.5
qy=0.5
nofp=100

printf -v num "%05i" $1
  
echo ./havoc --section=$num.hc --px=$px --py=$py --qx=$qx --qy=$qy --nofp=$nofp --var=$var 
./havoc --section=$num.hc --px=$px --py=$py --qx=$qx --qy=$qy --nofp=$nofp --var=$var 


