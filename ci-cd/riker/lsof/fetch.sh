#!/bin/bash

TOP="$(git rev-parse --show-toplevel)"
input_dir="${TOP}/ci-cd/inputs/scripts/lsof"

mkdir -p "$input_dir/dev"

git clone https://github.com/lsof-org/lsof.git "$input_dir/dev"
git -C "$input_dir/dev" checkout 005e014e1abdadb2493d8b3ce87b37a2c0a2351d
export BENCHMARK_INPUT_FILE="$(realpath "$input_dir/dev")"
(cd "$input_dir/dev" && "$input_dir/dev/Configure" -n linux)

