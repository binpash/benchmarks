#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
scripts_dir="${eval_dir}/scripts"

for bench in "$scripts_dir"/*; do
    "$bench/clean.sh" "$@"
done

