#!/bin/bash

printf -v num "%05i" $1
sed "s/_NUM_/$num/g" multiplot.plt > temp
sed "s/_NUM_/$num/g" temp > temp2.plt

if [ $# -eq 2 ];then
	if [ $2 = "g" ];then

	sed "s/_CASE1_/""/g" temp2.plt > temp3.plt	
	sed "s/_CASE2_/"#"/g" temp3.plt > temp4.plt

	else

	sed "s/_CASE2_/""/g" temp2.plt > temp3.plt	
	sed "s/_CASE1_/"#"/g" temp3.plt > temp4.plt

	fi
else

sed "s/_CASE2_/""/g" temp2.plt > temp3.plt	
sed "s/_CASE1_/"#"/g" temp3.plt > temp4.plt

fi


gnuplot "temp4.plt"

rm temp temp2.plt temp3.plt temp4.plt
