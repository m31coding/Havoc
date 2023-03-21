#!/bin/bash

if [ $# -ne 3 ]; then
echo "ERROR: three arguments needed."
exit 0
fi

j=0
for((i=$1;i<=$2;i+=$3))
do
echo "$i / $2"
printf -v num "%05i" $i
printf -v time "%.3f" $t
sed "s/_NUM_/$num/g" multiplot.plt > temp
j=$(($j+1)) 
gnuplot "temp"
done

rm temp 

ffmpeg -r 25 -i %05d.png -q 1 "multiplot_movie.mp4"

#rm *.png
