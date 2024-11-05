#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
scripts_dir="${eval_dir}/scripts"
input_dir="${eval_dir}/input"

test_db="$scripts_dir/sqlite/data/hello-riker.sqlite3"

"$input_dir/scripts/sqlite/dev/sqlite3" "$test_db" 'SELECT * FROM hello' | diff -q - <(echo riker)
echo riker/sqlite $?
