#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/covid-mts"
outputs_dir="${eval_dir}/outputs"
scripts_dir="${eval_dir}/scripts"
hashes_dir="${eval_dir}/hashes"

suffix=""
if [[ "$@" == *"--small"* ]]; then
    suffix="_small"
fi

if [[ "$@" == *"--generate"* ]]; then
    # give relative paths to md5sum
    (cd "$outputs_dir"; md5sum "outputs$suffix"/* > "$hashes_dir/outputs$suffix.md5sum")
    exit 0
fi

# give relative paths to md5sum
(cd "$outputs_dir"; md5sum --check --quiet --status "$hashes_dir/outputs$suffix.md5sum")
echo covid-mts$suffix $?
