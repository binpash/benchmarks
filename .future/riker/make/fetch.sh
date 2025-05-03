#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts/make"

mkdir -p "$input_dir"

tar_path="$input_dir/make-4.3.tar.gz"
wget http://ftp.gnu.org/gnu/make/make-4.3.tar.gz -O "$tar_path"
tar xzf "$tar_path" -C "$input_dir"
mv "$input_dir/make-4.3" "$input_dir/dev"
rm "$tar_path"