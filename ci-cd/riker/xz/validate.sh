#!/bin/bash

TOP="$(git rev-parse --show-toplevel)"
input_dir="${TOP}/ci-cd/inputs"

xz_bin="$input_dir/scripts/xz/dev/xz"

echo "Hello, Riker!" \
    | "$xz_bin" \
    | "$xz_bin" -d \
    | diff -q - <(echo Hello, Riker!)
echo riker/xz $?
