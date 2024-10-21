#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/log-analysis"
input_dir="${eval_dir}/input"
scripts_dir="${eval_dir}/scripts"
hashes_dir="${eval_dir}/hashes"
results_dir="${eval_dir}/results"
mkdir -p $results_dir

suffix=".full"
if [[ "$@" == *"--small"* ]]; then
    suffix=".small"
fi

cd $results_dir # md5sum computes paths relative to cd

if [[ "$@" == *"--generate"* ]]; then
    md5sum pcaps$suffix/* > $hashes_dir/pcaps$suffix.md5sum
    md5sum nginx$suffix/* > $hashes_dir/nginx$suffix.md5sum
fi

okay=0
if ! md5sum --check --quiet $hashes_dir/pcaps$suffix.md5sum; then
    okay=1
    echo "img_convert $suffix failed verification"
fi
if ! md5sum --check --quiet $hashes_dir/nginx$suffix.md5sum; then
    okay=1
    echo "to_mp3 $suffix failed verification"
fi

exit $okay
