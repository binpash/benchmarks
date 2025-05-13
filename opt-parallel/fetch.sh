#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
URL='https://atlas.cs.brown.edu/data'
eval_dir="${REPO_TOP}/opt-parallel"
input_dir="${eval_dir}/inputs"

size=full
for arg in "$@"; do
    case "$arg" in
        --small) size=small ;;
        --min)   size=min ;;
    esac
done
# for args
if [ "$size" = "min" ]; then
    if [ -d "$input_dir/ChessData_min" ]; then
        echo "Data already downloaded and extracted."
        exit 0
    fi
    mkdir -p "$input_dir/ChessData_min"
    cp min_inputs/* "$input_dir/ChessData_min"
    exit 0
fi

if [ "$size" = "small" ]; then
    if [ -d "$input_dir/ChessData_small" ]; then
        echo "Data already downloaded and extracted."
        exit 0
    fi
    mkdir -p "$input_dir/ChessData_small"
    unzip -q "${URL}/ChessData_small.zip" -d "$input_dir/ChessData_small"
    exit 0
fi


if [ -d "$input_dir/ChessData" ]; then
    echo "Data already downloaded and extracted."
    exit 0
fi
mkdir -p "$input_dir/ChessData_small"
unzip -q "${URL}/ChessData_small.zip" -d "$input_dir/ChessData_small"
exit 0