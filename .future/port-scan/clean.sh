#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
outputs_dir="${REPO_TOP}/port-scan/outputs"
input_dir="${REPO_TOP}/port-scan/inputs"

rm -rf "$outputs_dir"
rm -rf "$input_dir"
