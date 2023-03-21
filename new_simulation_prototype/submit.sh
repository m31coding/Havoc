#!/bin/bash

SIM_NAME=$(basename $PWD)

sed "s/SIM_NAME/$SIM_NAME/g" submit_prototype.sh > submit_final.sh

qsub submit_final.sh
