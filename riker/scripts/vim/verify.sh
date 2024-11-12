#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
scripts_dir="${eval_dir}/scripts"
input_dir="${eval_dir}/input"

# Call compiled binary to write an empty file with a randomly chosen path.
# Must use -u /dev/null to specify a blank config because there might not be a ~/.vimrc, which vim would complain about.
# Assert that the file was created.

vim="$input_dir/scripts/vim/dev/src/vim"
canary="$(mktemp --dry-run --tmpdir="$input_dir/scripts/vim")"

test ! -f "$canary"
echo riker/vim/canary-did-not-exist $?

"$vim" -u /dev/null -e -c "w $canary" -c quit

test -f "$canary"
echo riker/vim/created-canary-file $?
