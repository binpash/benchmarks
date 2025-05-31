#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/ci-cd/riker"
input_dir="${REPO_TOP}/ci-cd/inputs/scripts"

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_SCRIPT="$(realpath "$eval_dir/redis/build.sh")"
export BENCHMARK_INPUT_FILE="$(realpath "$input_dir/redis/dev/src")"
(cd "$input_dir/redis/dev/src" && $KOALA_SHELL "$eval_dir/redis/build.sh")


