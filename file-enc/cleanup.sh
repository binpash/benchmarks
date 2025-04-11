#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/file-enc"
results_dir="${eval_dir}/results"
input_dir="${eval_dir}/input"

rm -rf $results_dir
rm -rf $input_dir