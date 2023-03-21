#!/bin/bash

param1="MESH_OMEGA"
param2="MESH_PSI"

range1='0.0'
range2='0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1'

for k in $range1; do
for l in $range2; do

~/scratch/sim_voronoi/ini_new_sim

grep -v $param1 new_simulation/constants > new_simulation/temp
echo "$param1	$k" >> new_simulation/temp

grep -v $param2 new_simulation/temp > new_simulation/temp2
echo "$param2	$l" >> new_simulation/temp2

rm new_simulation/temp 
mv new_simulation/temp2 new_simulation/constants


name=$k"_"$l

rm -r $name
mv new_simulation $name 

cd $name
./submit.sh
cd ..

done
done


