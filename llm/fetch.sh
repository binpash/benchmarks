#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/llm"
scripts_dir="${eval_dir}/scripts"


for bench in "$scripts_dir"/*; do
    "$bench/fetch.sh" "$@"
done
