#!/bin/bash

TOP=$(git rev-parse --show-toplevel)

eval_dir="${TOP}/weather"
outputs_dir="${eval_dir}/outputs"
scripts_dir="${eval_dir}/scripts"
input_dir="${eval_dir}/inputs"

export LC_ALL=C

size=full
for arg in "$@"; do
    case "$arg" in
        --small) size="small" ;;
        --min)   size="min" ;;
    esac
done

KOALA_SHELL=${KOALA_SHELL:-bash}

export BENCHMARK_CATEGORY="weather"

echo "max-temp"
export input_file="${input_dir}/temperatures.$size.txt"
export statistics_dir="$outputs_dir/statistics.$size"

mkdir -p "$statistics_dir"

BENCHMARK_INPUT_FILE="$(realpath "$input_file")"
export BENCHMARK_INPUT_FILE

BENCHMARK_SCRIPT="$(realpath "${scripts_dir}/temp-analytics.sh")"
export BENCHMARK_SCRIPT

$KOALA_SHELL "$scripts_dir/temp-analytics.sh"

echo "$?"

echo "tuft-weather"
export BENCHMARK_SCRIPT="$scripts_dir/tuft-weather.sh"
export BENCHMARK_INPUT_FILE="$input_dir/tuft_weather.${size}.txt"

mkdir -p "$outputs_dir/$size"

$KOALA_SHELL "$BENCHMARK_SCRIPT" "$BENCHMARK_INPUT_FILE" "$size" > "$outputs_dir/$size/turf_weather.log"
echo "$?"

rm -rf "$outputs_dir/$size/plots" || true
mkdir -p "$outputs_dir/$size/plots"

if [ -d "$eval_dir/plots" ]; then
    mv "$eval_dir/plots"/* "$outputs_dir/$size/plots/"
fi

