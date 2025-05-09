#!/bin/bash

export PASH_SPEC_TOP=${PASH_SPEC_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

output_dir="$PASH_SPEC_TOP/report/output/fully_seq"
result_dir="$PASH_SPEC_TOP/results/fully_seq"

./run --target sh-only
mkdir -p $result_dir/base
mv $output_dir/* $result_dir/base


for window in 0 10 20 40
do
    ./run --target hs-only --window $window
    mkdir -p $result_dir/$window
    mv $output_dir/* $result_dir/$window
done
