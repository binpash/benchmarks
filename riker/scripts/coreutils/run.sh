#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts"
scripts_dir="${eval_dir}/scripts"

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
export BENCHMARK_SCRIPT="$(realpath "$scripts_dir/coreutils/build.sh")"
(cd "$input_dir/coreutils/dev" && $BENCHMARK_SHELL "$scripts_dir/coreutils/build.sh")