#!/bin/bash

#one parameter: view plot via gnuplot
#3 parameters: make movie: first file, last file, step

if [ $# -eq 0 ];then
echo "ERROR: one, or three arguments needed for particles.sh"
fi

if [ $# -eq 2 ];then
echo "ERROR: one, or three arguments needed for particles.sh"
fi


#one argument (view plot via gnuplot, e.g. particles.sh 20
if [ $# -eq 1 ];then

printf -v num "%05i" $1

sed "s/_NUM_/$num/g" particles.plt > temp.plt
sed "s/_CASE1_/#/g" temp.plt > temp2.plt
sed "s/_CASE2_//g" temp2.plt > temp3.plt

gnuplot "temp3.plt"

rm temp.plt temp2.plt temp3.plt

fi


#three arguments: make a movie, e.g: particles.sh 20 60 2 makes a movie of the files 00020.data until 00060.data with every second file
if [ $# -eq 3 ];then

j=0

for((i=$1;i<=$2;i+=$3))
do
echo "$i / $2"
printf -v num "%05i" $i
printf -v picNum "%05i" $j


sed "s/_NUM_/$num/g" particles.plt > temp.plt
sed "s/_CASE1_//g" temp.plt > temp2.plt
sed "s/_CASE2_/#/g" temp2.plt > temp3.plt
sed "s/_NUMBER_/$picNum/g" temp3.plt > temp4.plt

j=$(($j+1))

gnuplot "temp4.plt"

done


rm temp.plt temp2.plt temp3.plt temp4.plt

ffmpeg -r 5 -sameq -i %05dparticles.png "particles_movie.mp4"

rm *.png

fi






