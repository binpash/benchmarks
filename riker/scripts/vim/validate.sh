#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input"

# Call compiled binary to write an empty file with a randomly chosen path.
# Must use -u /dev/null to specify a blank config because there might not be a ~/.vimrc, which vim would complain about.
# Assert that the file was created.

vim="$input_dir/scripts/vim/dev/src/vim"
canary="$(mktemp --dry-run --tmpdir="$input_dir/scripts/vim")"

test ! -f "$canary" > /dev/null 2>&1
echo riker/vim/canary-did-not-exist $?

"$vim" -u /dev/null -e -c "w $canary" -c quit > /dev/null 2>&1

test -f "$canary" 2>&1
echo riker/vim/created-canary-file $?