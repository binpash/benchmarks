#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
scripts_dir="${eval_dir}/scripts"
input_dir="${eval_dir}/input"

lsof_bin="$input_dir/scripts/lsof/dev/lsof"

tmp_file=$(mktemp)
tail -f "$tmp_file" > /dev/null &
tail_pid=$!

sleep 1

"$lsof_bin" -p "$tail_pid" | grep "$tmp_file"

status=$?
kill "$tail_pid"
rm "$tmp_file"

echo riker/lsof $status
exit $status
