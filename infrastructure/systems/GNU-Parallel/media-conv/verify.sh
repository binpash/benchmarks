#!/bin/bash

hash_audio_dir() {
    local src_dir=$1
    for src in $src_dir/*; do
        got_hash=$(ffmpeg -i "$src" -map 0:a -f md5 - 2>/dev/null)
        echo $got_hash $(realpath "--relative-to=$src_dir" "$src")
    done
}

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/media-conv"
results_dir="${eval_dir}/results"
hashes_dir="${eval_dir}/hashes"

suffix=".full"
if [[ "$@" == *"--small"* ]]; then
    suffix=".small"
fi

if [[ "$@" == *"--generate"* ]]; then
    bench=to_mp3$suffix
    hash_audio_dir "$results_dir/$bench" > "$hashes_dir/$bench.md5sum"

    cd $results_dir
    bench=img_convert$suffix
    md5sum $bench/* > "$hashes_dir/$bench.md5sum"

    exit 0
fi


bench=to_mp3$suffix
hash_audio_dir "$results_dir/$bench" | diff -q "$hashes_dir/$bench.md5sum" -
echo $bench $?

cd $results_dir
bench=img_convert$suffix
md5sum --check --quiet --status $hashes_dir/$bench.md5sum
echo $bench $?

