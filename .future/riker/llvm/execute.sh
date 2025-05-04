#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts"
scripts_dir="${eval_dir}/scripts"

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_SCRIPT="$(realpath "$scripts_dir/llvm/build.sh")"
(cd "$input_dir/llvm/dev" && $KOALA_SHELL "$scripts_dir/llvm/build.sh")


