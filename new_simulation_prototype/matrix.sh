#!/bin/bash

if [ $# -eq 0 ];then
echo "ERROR: one, two or three arguments needed"
fi

if [ $# -eq 3 ]; then
for((i=$1;i<=$2;i+=$3))
do
echo "$i / $2"
printf -v num "%05i" $i
sed "s/_NUM_/$num/g" matrix.plt > matrix_temp.plt
gnuplot matrix_temp.plt
epstopdf $num.eps 
done
fi


