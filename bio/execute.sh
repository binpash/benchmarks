#!/bin/bash

# create bam files with regions
################### 1KG SAMPLES
IN="inputs/bio-full"
IN_NAME="inputs/bio-full/input.txt"
OUT="outputs"

for arg in "$@"; do
    case "$arg" in
        --small)
            IN_NAME="inputs/bio-small/input_small.txt" 
            IN="inputs/bio-small"
            ;;
        --min)   
            IN_NAME="inputs/bio-min/input_min.txt" 
            IN="inputs/bio-min"
            ;;
    esac
done

size=full
subset=false
for arg in "$@"; do
    case "$arg" in
    --small) size=full # small uses a subset of full inputs
    subset=true
    ;;
    --min) size=min ;;
    esac
done

# export SIZE="$size" # for PARAMS.sh
export SIZE=full

KOALA_SHELL="${KOALA_SHELL:-bash}"
export BENCHMARK_CATEGORY="bio"
export KOALA_SHELL

script_file="./scripts/bio.sh"
BENCHMARK_SCRIPT="$(realpath "$script_file")"
export BENCHMARK_SCRIPT

BENCHMARK_INPUT_FILE="$(realpath "$IN")"
export BENCHMARK_INPUT_FILE

$KOALA_SHELL "$script_file" "$IN" "$IN_NAME" "$OUT"

# Note: The 'data.sh' script must be run first
teraseq_script_names="data
run_dRNASeq
run_5TERA"

if [[ "$size" == "min" ]]; then
    teraseq_script_names="data"
fi

if [[ "$subset" == true ]]; then
teraseq_script_names="data
run_dRNASeq"
fi
BENCHMARK_INPUT_FILE="$(realpath "inputs/full")"
export BENCHMARK_INPUT_FILE
while IFS= read -r script; do
    script_file="./scripts/$script.sh"
    BENCHMARK_SCRIPT="$(realpath "$script_file")"
    export BENCHMARK_SCRIPT

    echo "$script"
    $KOALA_SHELL "$script_file"
    echo "$?"
done <<< "$teraseq_script_names"
