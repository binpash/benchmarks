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
outputs_dir="${eval_dir}/outputs"
hashes_dir="${eval_dir}/hashes"

suffix=".full"
for arg in "$@"; do
    case "$arg" in
        --small) suffix=".small" ;;
        --min) suffix=".min" ;;
    esac
done

if [[ " $* " == *" --generate "* ]]; then
    bench=to_mp3$suffix
    hash_audio_dir "$outputs_dir/$bench" > "$hashes_dir/$bench.md5sum"

    cd $outputs_dir
    bench=img_convert$suffix
    md5sum $bench/* > "$hashes_dir/$bench.md5sum"

    exit 0
fi


bench=to_mp3$suffix
hash_audio_dir "$outputs_dir/$bench" | diff -q "$hashes_dir/$bench.md5sum" -
echo $bench $?

cd $outputs_dir
bench=img_convert$suffix
md5sum --check --quiet --status $hashes_dir/$bench.md5sum
echo $bench $?