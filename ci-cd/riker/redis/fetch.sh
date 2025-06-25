#!/bin/bash

TOP="$(git rev-parse --show-toplevel)"
input_dir="${TOP}/ci-cd/inputs/scripts/redis"

mkdir -p "$input_dir/dev"

git clone https://github.com/redis/redis "$input_dir/dev"
git -C "$input_dir/dev" checkout d96f47cf06b1cc24b82109e0e87ac5428517525a
(cd "$input_dir/dev" && make .make-prerequisites)

