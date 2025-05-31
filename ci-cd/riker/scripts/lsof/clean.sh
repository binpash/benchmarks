#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/ci-cd/riker"
input_dir="${eval_dir}/inputs/scripts/lsof"

rm -rf "$input_dir/dev"
