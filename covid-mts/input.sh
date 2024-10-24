#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/covid-mts"
input_dir="${eval_dir}/input"

mkdir -p "$input_dir"

curl --insecure 'https://atlas-group.cs.brown.edu/data/covid-mts/in.csv.gz' | gunzip > "$input_dir/in.csv"

curl --insecure 'https://atlas-group.cs.brown.edu/data/covid-mts/in_small.csv.gz' | gunzip > "$input_dir/in_small.csv"
