#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/ci-cd/riker"
scripts_dir="${eval_dir}/scripts"
input_dir="${eval_dir}/inputs"

"$input_dir/scripts/memcached/dev/memcached" --version | diff -q - <(echo memcached 1.6.9)
echo riker/memcached $?
