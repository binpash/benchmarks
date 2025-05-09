#!/bin/bash

export PATH=$PATH:$HOME/.local/bin
export PASH_SPEC_TOP=${PASH_SPEC_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}
export PASH_TOP=${PASH_TOP:-$PASH_SPEC_TOP/deps/pash}

output_dir="$PASH_SPEC_TOP/report/output/dgsh"
download_dir="$PASH_SPEC_TOP/report/resources/dgsh"
result_dir="$PASH_SPEC_TOP/results/dgsh"

./run --target sh-only
rm -rf /tmp/*
mkdir -p $result_dir/base
mv $output_dir/* $result_dir/base

for window in 0 5 10 20 40
do
    ./run --target hs-only --window $window
    rm -rf /tmp/*
    mkdir -p $output_dir/$result_dir
    mv $output_dir/* $result_dir/$window
done
