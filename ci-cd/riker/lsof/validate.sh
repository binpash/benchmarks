#!/bin/bash

TOP="$(git rev-parse --show-toplevel)"
eval_dir="${TOP}/ci-cd/riker"
input_dir="${eval_dir}/inputs"

lsof_bin="$input_dir/scripts/lsof/dev/lsof"

tmp_file=$(mktemp)
tail -f "$tmp_file" > /dev/null &
tail_pid=$!
sleep 1

"$lsof_bin" -p "$tail_pid" 2>/dev/null | grep "$tmp_file" > /dev/null 2>&1
status=$?
kill "$tail_pid" >/dev/null 2>&1
rm "$tmp_file" >/dev/null 2>&1

echo riker/lsof $status
exit $status
