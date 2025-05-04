#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts"
scripts_dir="${eval_dir}/scripts"

KOALA_SHELL=${KOALA_SHELL:-bash}
export BENCHMARK_SCRIPT="$(realpath "$scripts_dir/xz/build.sh")"
(cd "$input_dir/xz/dev" && $KOALA_SHELL "$scripts_dir/xz/build.sh")


