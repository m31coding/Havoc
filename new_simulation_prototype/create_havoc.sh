#!/bin/bash
dir=$(pwd)
cd ../../code
make
cd $dir
cp ../../code/havoc ./
