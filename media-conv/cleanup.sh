#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
results_dir="${REPO_TOP}/media-conv/results"
input_dir="${REPO_TOP}/media-conv/input"

rm -rf $results_dir
rm -rf $input_dir

