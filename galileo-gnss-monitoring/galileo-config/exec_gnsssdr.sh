#!/bin/bash

if [ $2 -eq 0 ]
then
    cd ../build
    make -j
    cd ../install
fi

#/home/harshad/projects/gsoc_21_gnss-sdr/install/gnss-sdr --config-file=$1
#/home/maryam/projects/personal/gnss-sdr-modified/install/gnss-sdr --config-file=$1
/home/maryam/prefix/src/gnss-sdr/install/gnss-sdr --config-file=$1
