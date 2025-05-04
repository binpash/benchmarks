#!/bin/bash

# create bam files with regions
################### 1KG SAMPLES
IN="inputs"
IN_NAME="input.txt"
OUT="outputs"

for arg in "$@"; do
    case "$arg" in
        --small) IN_NAME="input_small.txt" ;;
        --min)   IN_NAME="input_min.txt" ;;
    esac
done

KOALA_SHELL="${KOALA_SHELL:-bash}"
export BENCHMARK_CATEGORY="bio"
export KOALA_SHELL

script_file="./scripts/bio.sh"
BENCHMARK_SCRIPT="$(realpath "$script_file")"
export BENCHMARK_SCRIPT

BENCHMARK_INPUT_FILE="$(realpath "$IN_NAME")"
export BENCHMARK_INPUT_FILE

$KOALA_SHELL "$script_file" "$IN" "$IN_NAME" "$OUT"
