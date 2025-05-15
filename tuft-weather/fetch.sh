#!/bin/bash

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/tuft-weather"
inputs_dir="$eval_dir/inputs"

mkdir -p "$inputs_dir"

size=full
for arg in "$@"; do
    case "$arg" in
        --small) size=small ;;
        --min)   size=min ;;
    esac
done

if [ "$size" = "full" ]; then
    python3 $eval_dir/scripts/generate_input.py $inputs_dir/inputs_$size.txt --size $size
elif [ "$size" = "small" ]; then
    python3 $eval_dir/scripts/generate_input.py $inputs_dir/inputs_$size.txt --size $size
else
    python3 $eval_dir/scripts/generate_input.py $inputs_dir/inputs_$size.txt --size $size
fi