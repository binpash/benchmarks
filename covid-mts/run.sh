#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/covid-mts"
input_dir="${eval_dir}/input"
outputs_dir="${eval_dir}/outputs"
scripts_dir="${eval_dir}/scripts"

suffix=""
if [[ "$@" == *"--small"* ]]; then
    suffix="_small"
fi

input_file="$input_dir/in$suffix.csv"
output_scoped="$outputs_dir/outputs$suffix"
mkdir -p "$output_scoped"

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
export BENCHMARK_CATEGORY="covid-mts"
export BENCHMARK_INPUT_FILE="$(realpath "$input_file")"

export BENCHMARK_SCRIPT="$(realpath "$scripts_dir/1.sh")"
$BENCHMARK_SHELL "$scripts_dir/1.sh" "$input_file" > "$output_scoped/1.out"
export BENCHMARK_SCRIPT="$(realpath "$scripts_dir/2.sh")"
$BENCHMARK_SHELL "$scripts_dir/2.sh" "$input_file" > "$output_scoped/2.out"
export BENCHMARK_SCRIPT="$(realpath "$scripts_dir/3.sh")"
$BENCHMARK_SHELL "$scripts_dir/3.sh" "$input_file" > "$output_scoped/3.out"
export BENCHMARK_SCRIPT="$(realpath "$scripts_dir/4.sh")"
$BENCHMARK_SHELL "$scripts_dir/4.sh" "$input_file" > "$output_scoped/4.out"
