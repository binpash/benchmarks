#!/bin/bash

TOP="$(git rev-parse --show-toplevel)"
eval_dir="${TOP}/ci-cd/riker"
input_dir="${TOP}/ci-cd/inputs/scripts"

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_SCRIPT="$(realpath "$eval_dir/vim/build.sh")"
export BENCHMARK_INPUT_FILE="$(realpath "$input_dir/vim/dev")"
(cd "$input_dir/vim/dev" && $KOALA_SHELL "$eval_dir/vim/build.sh")


