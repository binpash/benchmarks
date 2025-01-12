#!/bin/bash

# create bam files with regions
################### 1KG SAMPLES
IN=inputs
IN_NAME=input.txt
OUT=outputs

if [[ "$@" == *"--small"* ]]; then
    IN_NAME=input_small.txt
fi

export BENCHMARK_CATEGORY="bio"
BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}

script_file=./scripts/bio.sh 
export BENCHMARK_SCRIPT="$(realpath "$script_file")"
export BENCHMARK_INPUT_FILE="$(realpath "$IN_NAME")"

time $BENCHMARK_SHELL "$script_file" "$IN" "$IN_NAME" "$OUT"
