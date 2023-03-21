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
printf -v picNum "%05i" $j
sed "s/_GRAPH_/$num.color2D/g" movie.plt > temp
sed "s/_FACES_/$num.faces/g" temp > temp2.plt
sed "s/_NUMBER_/$picNum/g" temp2.plt > temp3.plt
j=$(($j+1)) 
gnuplot "temp3.plt"
done

rm temp temp2.plt temp3.plt

ffmpeg -r 25 -i %05d.png -q 1 "color_movie.mp4"

#rm *.png
