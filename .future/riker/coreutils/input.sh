#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts/coreutils"

mkdir -p "$input_dir"

tar_path="$input_dir/coreutils-8.32.tar.gz"
wget https://ftp.gnu.org/gnu/coreutils/coreutils-8.32.tar.gz -O "$tar_path"
tar xzf "$tar_path" -C "$input_dir"
mv "$input_dir/coreutils-8.32" "$input_dir/dev"
rm "$tar_path"