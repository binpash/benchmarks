#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/covid-mts"
input_dir="${eval_dir}/inputs"

mkdir -p "$input_dir"

for arg in "$@"; do
    if [ "$arg" = "--min" ]; then
        cp min_inputs/in_min.csv "$input_dir/in_min.csv"
        exit 0
    fi
    if [ "$arg" = "--small" ]; then
        if [ ! -f "$input_dir/in_small.csv" ]; then
            curl --insecure 'https://atlas-group.cs.brown.edu/data/covid-mts/in_small.csv.gz' | gunzip > "$input_dir/in_small.csv"
            exit 0
        else
            exit 0
        fi
    fi
done

if [ -f "$input_dir/in.csv" ]; then
    exit 0
fi
curl --insecure 'https://atlas-group.cs.brown.edu/data/covid-mts/in_full.csv.gz' | gunzip > "$input_dir/in.csv"
