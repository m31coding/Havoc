#!/bin/bash

if [ $# -eq 0 ];then
echo "ERROR: one, two or three arguments needed for convert.sh"
fi

#one argument: convert the file. e.g convert.sh 3 
if [ $# -eq 1 ]; then
printf -v num "%05i" $1
./havoc --convert="$num.hc" --hybrid --color_only
fi

#two arguments: convert all files between the arguments. e.g convert.sh 3 10
if [ $# -eq 2 ]; then
for((i=$1;i<=$2;i++))
do
echo "$i / $2"
printf -v num "%05i" $i
./havoc --convert="$num.hc" --hybrid --color_only
done
fi

#three arguments: convert all files between the arguments
#but convert every $3.th file only. e.g convert 4 10 2
if [ $# -eq 3 ]; then
for((i=$1;i<=$2;i+=$3))
do
echo "$i / $2"
printf -v num "%05i" $i
./havoc --convert="$num.hc" --hybrid --color_only
done
fi


