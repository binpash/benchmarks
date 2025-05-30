#!/bin/bash

KOALA_SHELL="${KOALA_SHELL:-bash}"

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/tuft-weather"
scripts_dir="$eval_dir/scripts"
inputs_dir="$eval_dir/inputs"
outputs_dir="$eval_dir/outputs"

rm -rf "$outputs_dir" || true
mkdir -p "$outputs_dir"

size="full"
for arg in "$@"; do
    case "$arg" in
        --small) size="small" ;;
        --min)   size="min" ;;
    esac
done

export BENCHMARK_CATEGORY="tuft-weather"
export BENCHMARK_SCRIPT="$scripts_dir/tuft-weather.sh"
export BENCHMARK_INPUT_FILE="$inputs_dir/inputs_${size}.txt"

mkdir -p "$outputs_dir/$size"

$KOALA_SHELL "$BENCHMARK_SCRIPT" "$BENCHMARK_INPUT_FILE" "$size" > "$outputs_dir/$size/turf_weather.log"
echo "$?"

rm -rf "$outputs_dir/$size/plots" || true
mkdir -p "$outputs_dir/$size/plots"

if [ -d "$eval_dir/plots" ]; then
    mv "$eval_dir/plots"/* "$outputs_dir/$size/plots/"
fi