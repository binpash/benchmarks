#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts/memcached"

mkdir -p "$input_dir/dev"

git clone https://github.com/memcached/memcached.git "$input_dir/dev"
git -C "$input_dir/dev" checkout c472369fed5981ba8c004d426cee62d5165c47ca 

(cd "$input_dir/dev" && "$input_dir/dev/autogen.sh")

(cd "$input_dir/dev" && "$input_dir/dev/configure" --disable-dependency-tracking)


