#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts/llvm"

mkdir -p "$input_dir/dev"

git clone https://github.com/llvm/llvm-project.git "$input_dir/dev"
git -C "$input_dir/dev" checkout d28af7c654d8db0b68c175db5ce212d74fb5e9bc

