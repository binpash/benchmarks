#!/bin/bash

TOP="$(git rev-parse --show-toplevel)"
input_dir="${TOP}/ci-cd/inputs"

"$input_dir/scripts/memcached/dev/memcached" --version | diff -q - <(echo memcached 1.6.9)
echo riker/memcached $?
