#!/bin/bash

TOP="$(git rev-parse --show-toplevel)"
eval_dir="${TOP}/ci-cd/riker"
input_dir="${eval_dir}/inputs"

"$input_dir/scripts/lua/dev/src/lua" -e 'print("Hello, Riker!")' | diff -q - <(echo "Hello, Riker!")
echo riker/lua $?
