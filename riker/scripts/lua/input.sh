#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts/lua"

mkdir -p "$input_dir/dev"

wget https://www.lua.org/ftp/lua-5.4.3.tar.gz -O "$input_dir"
tar xzf "$input_dir/lua-5.4.3.tar.gz" -C "$input_dir/dev"
