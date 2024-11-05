#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
scripts_dir="${eval_dir}/scripts"
input_dir="${eval_dir}/input"

vim="$input_dir/scripts/vim/dev/src/vim"
canary=$(mktemp --dry-run --tmpdir="$input_dir/scripts/vim")

"$vim" -u /dev/null -e -c "w $canary" -c quit
test -f $canary
echo riker/vim $?
