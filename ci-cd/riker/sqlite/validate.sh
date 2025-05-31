#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/ci-cd/riker"
input_dir="${REPO_TOP}/ci-cd/inputs"

test_db="$eval_dir/sqlite/data/hello-riker.sqlite3"

"$input_dir/scripts/sqlite/dev/sqlite3" "$test_db" 'SELECT * FROM hello' | diff -q - <(echo riker)
echo riker/sqlite $?
