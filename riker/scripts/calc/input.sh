#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts/calc"

mkdir -p "$input_dir/dev"

git clone https://github.com/lcn2/calc.git "$input_dir/dev"
git -C "$input_dir/dev" checkout b0f19c10110feacd630b3b48fb4c1080e80bae28

