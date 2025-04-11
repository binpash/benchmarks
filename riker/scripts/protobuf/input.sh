#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts/protobuf"

mkdir -p "$input_dir/dev"

git clone https://github.com/protocolbuffers/protobuf.git "$input_dir/dev"
git -C "$input_dir/dev" checkout 909a0f36a10075c4b4bc70fdee2c7e32dd612a72

