#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/covid-mts"
input_dir="${eval_dir}/inputs"

mkdir -p "$input_dir"

for arg in "$@"; do
    if [ "$arg" = "--min" ]; then
        curl --insecure 'https://atlas-group.cs.brown.edu/data/covid-mts/in_min.csv.gz' | gunzip > "$input_dir/in_min.csv"
        exit 0
    fi
    if [ "$arg" = "--small" ]; then
        curl --insecure 'https://atlas-group.cs.brown.edu/data/covid-mts/in_small.csv.gz' | gunzip > "$input_dir/in_small.csv"
        exit 0
    fi
done

curl --insecure 'https://atlas-group.cs.brown.edu/data/covid-mts/in.csv.gz' | gunzip > "$input_dir/in.csv"
