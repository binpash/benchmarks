#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/ci-cd"

for bench in "$eval_dir"/*; do
    if [ ! -d "$bench" ]; then
        continue
    fi
    "$bench/install.sh" "$@"
done
