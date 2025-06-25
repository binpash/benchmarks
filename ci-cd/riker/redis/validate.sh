#!/bin/bash

TOP="$(git rev-parse --show-toplevel)"
input_dir="${TOP}/ci-cd/inputs"

"$input_dir/scripts/redis/dev/src/redis-cli" --version | diff - <(echo "redis-cli 255.255.255 (git:d96f47cf)")
echo riker/redis $?
