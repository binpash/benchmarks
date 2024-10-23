#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
input_dir="${REPO_TOP}/max-temp/input"
results_dir="${REPO_TOP}/max-temp/results"

rm -rf "$input_dir"
rm -rf "$results_dir"

