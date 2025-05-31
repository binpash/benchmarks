#!/bin/bash

TOP="$(git rev-parse --show-toplevel)"
eval_dir="${TOP}/ci-cd/riker"
input_dir="${eval_dir}/inputs"

"$input_dir/scripts/memcached/dev/memcached" --version | diff -q - <(echo memcached 1.6.9)
echo riker/memcached $?
