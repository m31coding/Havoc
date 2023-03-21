#!/bin/bash

if [ $# -ne 1 ];then
echo "ERROR in graph.sh: need a filenumber as argument, e.g. ./graph.sh 100"
exit
fi

cd graph
./clean
cd ..

printf -v num "%05i" $1

./havoc --graph=$num.hc

cd graph
gnuplot "graph_ghost.plt"




