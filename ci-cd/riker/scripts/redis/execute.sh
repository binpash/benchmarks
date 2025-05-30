#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/inputs/scripts"
scripts_dir="${eval_dir}/scripts"

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_SCRIPT="$(realpath "$scripts_dir/redis/build.sh")"
export BENCHMARK_INPUT_FILE="$(realpath "$input_dir/redis/dev/src")"
(cd "$input_dir/redis/dev/src" && $KOALA_SHELL "$scripts_dir/redis/build.sh")


