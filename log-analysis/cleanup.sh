#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
input_dir="${REPO_TOP}/log-analysis/inputs"
outputs_dir="${REPO_TOP}/log-analysis/outputs"

rm -rf "$input_dir"
rm -rf "$outputs_dir"
