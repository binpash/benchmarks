#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/file-enc"
results_dir="${eval_dir}/outputs"
input_dir="${eval_dir}/inputs"

rm -rf "$results_dir"
rm -rf "$input_dir"