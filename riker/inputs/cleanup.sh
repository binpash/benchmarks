#/bin/bash

cd ..
benchmarks="lsof memcached redis sqlite vim xz"

for benchmark in $benchmarks; do
    cd $benchmark
    git clean -fdx
    git reset --hard
    cd ..
done

cd lua
make clean