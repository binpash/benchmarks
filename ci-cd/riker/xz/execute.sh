#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/ci-cd/riker"
input_dir="${REPO_TOP}/ci-cd/inputs/scripts"

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_SCRIPT="$(realpath "$eval_dir/xz/build.sh")"
export BENCHMARK_INPUT_FILE="$(realpath "$input_dir/xz/dev")"
(cd "$input_dir/xz/dev" && $KOALA_SHELL "$eval_dir/xz/build.sh")


