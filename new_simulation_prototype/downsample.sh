#!/bin/bash

#Latex picture length: 6.2 inches
#printer: 300 dpi => 1860 pixel => ~ 2000 pixel 

if [ $# -ne 3 ]; then
echo "ERROR: three arguments needed."
exit 0
fi

if [ $# -eq 3 ]; then
for((i=$1;i<=$2;i+=$3))
do
echo "$i / $2"
printf -v num "%05i" $i
./havoc --downsampling=$num.hc --pixel_per_length=200 --var=density --px=-0.49 --py=-0.49 --qx=0.49 --qy=0.49
done
fi

