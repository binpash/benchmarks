#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
input_dir="${REPO_TOP}/ci-cd/inputs/scripts/sqlite"

mkdir -p "$input_dir/dev"

git clone https://github.com/sqlite/sqlite "$input_dir/dev"
git -C "$input_dir/dev" checkout c1cace0832fa2af5ab8315e217d708c09d586425
(cd "$input_dir/dev" && "$input_dir/dev/configure")

