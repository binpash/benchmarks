#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
outputs_dir="${REPO_TOP}/dpt/outputs"
input_dir="${REPO_TOP}/dpt/input"

rm -rf "$outputs_dir"
rm -rf "$input_dir"