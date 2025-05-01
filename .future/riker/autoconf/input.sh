#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts/autoconf"

mkdir -p "$input_dir"

tar_path="$input_dir/autoconf-2.69.tar.gz"
wget https://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz -O "$tar_path"
tar xzf "$tar_path" -C "$input_dir"
mv "$input_dir/autoconf-2.69" "$input_dir/dev"
rm "$tar_path"

