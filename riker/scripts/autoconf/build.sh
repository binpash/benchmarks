#!/bin/sh

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts/autoconf"
temp_dir="${eval_dir}/input/scripts/temp"

mkdir -p "$temp_dir"
./configure --prefix="$temp_dir"
make
make install