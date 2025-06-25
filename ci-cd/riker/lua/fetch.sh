#!/bin/bash

TOP="$(git rev-parse --show-toplevel)"
input_dir="${TOP}/ci-cd/inputs/scripts/lua"

mkdir -p "$input_dir"

tar_path="$input_dir/lua-5.4.3.tar.gz"
wget https://www.lua.org/ftp/lua-5.4.3.tar.gz -O "$tar_path"
tar xzf "$tar_path" -C "$input_dir"
mv "$input_dir/lua-5.4.3" "$input_dir/dev"
rm "$tar_path"