#!/bin/bash

OPTION="color_vorticity"

if [ $# -eq 0 ];then
echo "ERROR: one, two or three arguments needed for convert.sh"
fi

#one argument: color the file. e.g convert.sh 3 
if [ $# -eq 1 ]; then
printf -v num "%05i" $1
./havoc --$OPTION="$num.hc" 
fi

#two arguments: coler all files between the arguments. e.g convert.sh 3 10
if [ $# -eq 2 ]; then
for((i=$1;i<=$2;i++))
do
echo "$i / $2"
printf -v num "%05i" $i
./havoc --$OPTION="$num.hc"
done
fi

#three arguments: coler all files between the arguments
#but convert every $3.th file only. e.g convert 4 10 2
if [ $# -eq 3 ]; then
for((i=$1;i<=$2;i+=$3))
do
echo "$i / $2"
printf -v num "%05i" $i
./havoc --$OPTION="$num.hc" 
done
fi


