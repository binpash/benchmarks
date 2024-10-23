#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
results_dir="${REPO_TOP}/log-analysis/input"
results_dir="${REPO_TOP}/log-analysis/results"

rm -rf $input_dir
rm -rf $results_dir
