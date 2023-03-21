#!/bin/bash

if [ $# -eq 0 ]; then

rm *.hc *.data *.density2D *.pressure2D *.faces *.color2D *.png sim_info 
rm *.section *.sedov
rm Sim.*
rm *.ghosts *.ghost_faces
rm *.obstacle

fi

#two arguments: clean all files between the arguments. e.g clean.sh 3 10
if [ $# -eq 2 ]; then
for((i=$1;i<=$2;i++))
do
echo "$i / $2"
printf -v num "%05i" $i
rm $num.*
done
fi

