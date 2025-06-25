#!/bin/bash

TOP="$(git rev-parse --show-toplevel)"
eval_dir="${TOP}/ci-cd/riker"
input_dir="${TOP}/ci-cd/inputs/scripts"

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_SCRIPT="$(realpath "$eval_dir/lsof/build.sh")"
export BENCHMARK_INPUT_FILE="$(realpath "$input_dir/lsof/dev")"
(cd "$input_dir/lsof/dev" && $KOALA_SHELL "$eval_dir/lsof/build.sh")

