#!/bin/bash

if [ $# -ne 1 ]; then
echo ERROR in cut_convert.sh: one argument is needed, e.g. 00010.hc
exit
fi

xmin=-0.49
ymin=-0.49
xmax=0.49
ymax=0.49

printf -v num "%05i" $1

echo havoc --convert=$num.hc --cut --xmin=$xmin --ymin=$ymin --xmax=$xmax --ymax=$ymax  
./havoc --convert=$num.hc --cut --xmin=$xmin --ymin=$ymin --xmax=$xmax --ymax=$ymax  

