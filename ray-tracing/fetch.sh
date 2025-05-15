#!/bin/bash
set -e

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/ray-tracing"
input_dir="$eval_dir/inputs"
mkdir -p "$input_dir"

size="full"
for arg in "$@"; do
    case "$arg" in
    --small) size="small" ;;
    --min) size="min" ;;
    esac
done
input_dir="$input_dir/$size"
mkdir -p "$input_dir"
if [ "$size" = "small" ]; then
    N=400000
elif [ "$size" = "min" ]; then
    N=4000
else
    N=40000000
fi

python3 generate_inputs.py "$input_dir" "$N"
