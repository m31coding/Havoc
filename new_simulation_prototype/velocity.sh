#!/bin/bash

scale=0.02

printf -v num "%05i" $1
sed "s/_FILE_/"$num".data/g" velocity.plt > temp
sed "s/_SCALE_/$scale/g" temp > temp2.plt

gnuplot "temp2.plt"

rm temp temp2.plt

