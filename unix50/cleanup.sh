#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/unix50"
results_dir="${eval_dir}/results"

rm -rf $results_dir
