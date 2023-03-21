#!/bin/bash

DIR=$PWD
SRC_DIR="../source_code"
SRC_DIR_SIM="./src"

rm -rf $SRC_DIR_SIM
cp -r $SRC_DIR $SRC_DIR_SIM
cp ini.cpp $SRC_DIR_SIM
cd $SRC_DIR_SIM
make
cd $DIR
cp $SRC_DIR_SIM/havoc $DIR/havoc
