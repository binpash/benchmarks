#!/bin/bash

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/tuft-weather"
scripts_dir="$eval_dir/scripts"
inputs_dir="$eval_dir/inputs"
outputs_dir="$eval_dir/outputs"
mkdir -p "$outputs_dir"
source $eval_dir/venv/bin/activate
rm -r $eval_dir/plots
size=full
for arg in "$@"; do
    case "$arg" in
        --small) size=small ;;
        --min)   size=min ;;
    esac
done

export BENCHMARK_CATEGORY="tuft-weather"

export BENCHMARK_SCRIPT="$scripts_dir/tuft-weather.sh"
export BENCHMARK_INPUT_FILE="$inputs_dir/inputs_$size.txt"
mkdir -p "$outputs_dir/$size"
$BENCHMARK_SHELL "$scripts_dir/tuft-weather.sh" "$inputs_dir/inputs_$size.txt" $size > "$outputs_dir/$size/turf_weather.log"
echo "$?"

rm -rf $eval_dir/outputs/$size/plots
mkdir -p $eval_dir/outputs/$size/plots
mv $eval_dir/plots $eval_dir/outputs/$size/plots/