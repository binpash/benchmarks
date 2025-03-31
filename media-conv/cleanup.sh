#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
outputs_dir="${REPO_TOP}/media-conv/outputs"
input_dir="${REPO_TOP}/media-conv/inputs"

rm -rf "$outputs_dir"
rm -rf "$input_dir"

