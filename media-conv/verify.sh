#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/media-conv"
results_dir="${eval_dir}/results"
hashes_dir="${eval_dir}/hashes"

suffix=".full"
if [[ "$@" == *"--small"* ]]; then
    suffix=".small"
fi

cd $results_dir # md5sum computes paths relative to cd

if [[ "$@" == *"--generate"* ]]; then
    md5sum img_convert$suffix/* > $hashes_dir/img_convert$suffix.md5sum
    md5sum to_mp3$suffix/* > $hashes_dir/to_mp3$suffix.md5sum
    exit 0
fi

bench=img_convert$suffix
md5sum --check --quiet --status $hashes_dir/$bench.md5sum
echo $bench $?

bench=to_mp3$suffix
md5sum --check --quiet --status $hashes_dir/$bench.md5sum
echo $bench $?
