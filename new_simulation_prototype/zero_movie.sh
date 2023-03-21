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
#t=$(head -1 "$num.data" | awk '{print $2}')
printf -v time "%.3f" $t
./havoc --graph=$num.hc
cd graph
sed "s/_GRAPH_/$num.density2D/g" movie_yt.plt > temp
sed "s/_FACES_/$num.faces/g" temp > temp2.plt
sed "s/_NUMBER_/$picNum/g" temp2.plt > temp3.plt
sed "s/_TIME_/$time/g" temp3.plt > temp4.plt
j=$(($j+1)) 
gnuplot "temp4.plt"
mv $num.png ../
cd ..
done

rm temp temp2.plt temp3.plt temp4.plt

ffmpeg -r 25 -i %05d.png -q 1 "zero_movie.mp4"

#rm *.png
