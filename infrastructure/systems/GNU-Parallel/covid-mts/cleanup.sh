#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/covid-mts"
outputs_dir="${eval_dir}/outputs"
input_dir="${eval_dir}/input"

rm -rf "$outputs_dir"
rm -rf "$input_dir"
