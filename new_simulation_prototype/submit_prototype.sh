#!/bin/bash
DIR=$PWD

#$ -N "SIM_NAME"
#$ -M emailaddr
# -pe sol 8
#$ -cwd
#$ -q workstations
#$ -l h_vmem=1G

cd $DIR
hostname
date
./havoc 
date

cp 00000.hc ./start/

