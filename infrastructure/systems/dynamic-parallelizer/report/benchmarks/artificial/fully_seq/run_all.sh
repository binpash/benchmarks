#!/bin/bash

./setup.sh
# rm /mydata/input.txt
time ./run --target sh-only
mkdir -p /mydata/results/output/artificial/fully_seq/base
rm /mydata/input.txt
mv /mydata/dynamic-parallelizer/report/output/artificial/fully_seq/* /mydata/results/output/artificial/fully_seq/base

for window in 0 
do
    ./setup.sh
    time ./run --target hs-only --window $window
    rm /mydata/input.txt
    mkdir -p /mydata/results/output/artificial/fully_seq/$window
    mv /mydata/dynamic-parallelizer/report/output/artificial/fully_seq/* /mydata/results/output/artificial/fully_seq/$window
done

