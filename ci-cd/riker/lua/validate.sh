#!/bin/bash

TOP="$(git rev-parse --show-toplevel)"
input_dir="${TOP}/ci-cd/inputs"

"$input_dir/scripts/lua/dev/src/lua" -e 'print("Hello, Riker!")' | diff -q - <(echo "Hello, Riker!")
echo riker/lua $?
