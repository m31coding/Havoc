#!/bin/bash

printf -v num "%05i" $1

sed "s/_NUM_/$num/g" density2D.plt > temp.plt

if [ $# -eq 2 ];then
	if [ $2 = "g" ];then

	sed "s/_CASE1_/""/g" temp.plt > temp2.plt	
	sed "s/_CASE2_/"#"/g" temp2.plt > temp3.plt
	sed "s/_CASE3_/"#"/g" temp3.plt > temp4.plt

	else

	if [ $2 = "gg" ];then

	sed "s/_CASE1_/"#"/g" temp.plt > temp2.plt	
	sed "s/_CASE2_/"#"/g" temp2.plt > temp3.plt
	sed "s/_CASE3_/""/g" temp3.plt > temp4.plt

	else

	sed "s/_CASE1_/"#"/g" temp.plt > temp2.plt	
	sed "s/_CASE2_/""/g" temp2.plt > temp3.plt
	sed "s/_CASE3_/"#"/g" temp3.plt > temp4.plt
	
	fi
	fi

else
	sed "s/_CASE1_/"#"/g" temp.plt > temp2.plt	
	sed "s/_CASE2_/""/g" temp2.plt > temp3.plt
	sed "s/_CASE3_/"#"/g" temp3.plt > temp4.plt
fi


gnuplot "temp4.plt"

rm temp.plt temp2.plt temp3.plt temp4.plt

