#!/bin/bash
DIR=$PWD

#$ -N "Job"
#$ -M emailaddr
# -pe sol 8
#$ -cwd
#$ -q workstations

cd $DIR
hostname
date
./convert.sh 0 1000 1 
date


