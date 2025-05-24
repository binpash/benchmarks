#!/bin/bash

# create bam files with regions
################### 1KG SAMPLES
IN="inputs/full"
IN_NAME="inputs/full/input.txt"
OUT="outputs"

for arg in "$@"; do
    case "$arg" in
        --small)
            IN_NAME="inputs/small/input_small.txt" 
            IN="inputs/small"
            ;;
        --min)   
            IN_NAME="inputs/min/input_min.txt" 
            IN="inputs/min"
            ;;
    esac
done

KOALA_SHELL="${KOALA_SHELL:-bash}"
export BENCHMARK_CATEGORY="bio"
export KOALA_SHELL

script_file="./scripts/bio.sh"
BENCHMARK_SCRIPT="$(realpath "$script_file")"
export BENCHMARK_SCRIPT

BENCHMARK_INPUT_FILE="$(realpath "$IN")"
export BENCHMARK_INPUT_FILE

$KOALA_SHELL "$script_file" "$IN" "$IN_NAME" "$OUT"
