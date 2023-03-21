#!/bin/bash

#rename png files

n_shift=-100
num=0

for((i=$1;i<=$2;i++))
do
num=$(echo "$i + $n_shift" | bc)
printf -v num_old "%05d" $i
printf -v num_new "%05d" $num

mv "$num_old.png" "$num_new.png"

done


