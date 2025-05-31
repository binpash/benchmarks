#!/bin/bash

TOP="$(git rev-parse --show-toplevel)"
eval_dir="${TOP}/ci-cd/riker"
input_dir="${eval_dir}/inputs"

"$input_dir/scripts/redis/dev/src/redis-cli" --version | diff - <(echo "redis-cli 255.255.255 (git:d96f47cf)")
echo riker/redis $?
