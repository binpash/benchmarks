#!/bin/bash

cd ..
benchmarks="lsof lua memcached redis sqlite vim xz"

for benchmark in $benchmarks; do
    ./$benchmark"_build.sh"
done