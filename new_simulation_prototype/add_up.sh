#!/bin/bash

var="mass_std_dev_square"
#var="turbulence_kinetic_energy"
#var="kinetic_energy"
#var="abs_vorticity"
#var="energy"
#var="kinetic_energy_y"

rm $var\_data
touch $var\_data

if [ $# -eq 3 ]; then
for((i=$1;i<=$2;i+=$3))
do
echo "$i / $2"
printf -v num "%05i" $i
./havoc --add_up=$num.hc --$var >> $var\_data 
done
else
echo "Error: three arguments needed"
fi

