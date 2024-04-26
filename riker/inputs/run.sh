#/bin/bash

cd ..
benchmarks="lsof lua memcached redis sqlite vim xz"
tests="test"

for benchmark in $tests; do
    source $benchmark"_venv/bin/activate"
    ./$benchmark"_build.sh"
    deactivate
done