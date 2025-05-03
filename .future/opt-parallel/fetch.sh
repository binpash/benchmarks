#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/opt-parallel"
inputs_dir="${eval_dir}/inputs"

mkdir -p $inputs_dir
cd $inputs_dir || exit 1
# for args
for arg in "$@"; do
    if [[ "$arg" == "--small" ]]; then
        mkdir -p $inputs_dir/ChessData
        cp $eval_dir/temp_inputs/* $inputs_dir/ChessData/
        exit 0
    fi
    if [[ "$arg" == "--min" ]]; then
        mkdir -p $inputs_dir/ChessData
        cp $eval_dir/temp_inputs/* $inputs_dir/ChessData/
        exit 0
    fi
done
if [[ -d $inputs_dir/ChessData ]]; then
    echo "Data already downloaded and extracted."
    exit 0
fi
git clone https://github.com/rozim/ChessData.git

# TODO: Add small inputs and link for download