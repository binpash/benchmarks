#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/ci-cd/riker"
input_dir="${eval_dir}/inputs/scripts"
scripts_dir="${eval_dir}/scripts"

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_SCRIPT="$(realpath "$scripts_dir/lsof/build.sh")"
export BENCHMARK_INPUT_FILE="$(realpath "$input_dir/lsof/dev")"
(cd "$input_dir/lsof/dev" && $KOALA_SHELL "$scripts_dir/lsof/build.sh")

