#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts/sqlite"

mkdir -p "$input_dir/dev"

git clone https://github.com/sqlite/sqlite "$input_dir/dev"
git -C "$input_dir/dev" checkout c1cace0832fa2af5ab8315e217d708c09d586425
(cd "$input_dir/dev" && "$input_dir/dev/configure")

