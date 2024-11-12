#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
scripts_dir="${eval_dir}/scripts"
input_dir="${eval_dir}/input"

"$input_dir/scripts/redis/dev/src/redis-cli" --version | diff - <(echo "redis-cli 255.255.255 (git:d96f47cf)")
echo riker/redis $?
