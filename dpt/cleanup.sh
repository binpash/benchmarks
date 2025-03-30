#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
results_dir="${REPO_TOP}/dpt/results"
input_dir="${REPO_TOP}/dpt/input"

rm -rf "$results_dir"
rm -rf "$input_dir"
