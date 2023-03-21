#!/bin/bash

var="density"
scale=0.0001

printf -v num "%05i" $1
sed "s/_FILE_/"$num"_$var.gradients/g" gradients.plt > temp
sed "s/_SCALE_/$scale/g" temp > temp2.plt

gnuplot "temp2.plt"

rm temp temp2.plt

