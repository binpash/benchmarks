#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts"
scripts_dir="${eval_dir}/scripts"

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
(cd "$input_dir/redis/dev/src" && $BENCHMARK_SHELL "$scripts_dir/redis/build.sh")


