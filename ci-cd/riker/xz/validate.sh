#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/ci-cd/riker"
input_dir="${REPO_TOP}/ci-cd/inputs"

xz_bin="$input_dir/scripts/xz/dev/xz"

echo "Hello, Riker!" \
    | "$xz_bin" \
    | "$xz_bin" -d \
    | diff -q - <(echo Hello, Riker!)
echo riker/xz $?
