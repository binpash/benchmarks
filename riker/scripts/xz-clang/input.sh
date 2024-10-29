#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts/xz-clang"

mkdir -p "$input_dir/dev"

git clone https://github.com/xz-mirror/xz "$input_dir/dev"
git -C "$input_dir/dev" checkout 2327a461e1afce862c22269b80d3517801103c1b

