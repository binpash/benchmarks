#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
scripts_dir="${eval_dir}/scripts"
input_dir="${eval_dir}/input"

"$input_dir/scripts/lua/dev/src/lua" -e 'print("Hello, Riker!")' | diff -q - <(echo "Hello, Riker!")
echo riker/lua $?
