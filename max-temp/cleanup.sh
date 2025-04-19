#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
input_dir="${REPO_TOP}/max-temp/inputs"
outputs_dir="${REPO_TOP}/max-temp/outputs"

rm -rf "${input_dir}"
rm -rf "${outputs_dir}"
