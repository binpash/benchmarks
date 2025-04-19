#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
input_dir="${REPO_TOP}/tuft-weather/inputs"
outputs_dir="${REPO_TOP}/tuft-weather/outputs"

rm -rf "${input_dir}"
rm -rf "${outputs_dir}"
rm *.txt
rm *.png