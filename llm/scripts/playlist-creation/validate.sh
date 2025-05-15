#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/llm/scripts/playlist-creation"
hashes_dir="${eval_dir}/hashes"

suffix=""
generate=false
for arg in "$@"; do
    case "$arg" in
        --generate) generate=true ;;
        --small) suffix=".small" ;;
        --min) suffix=".min" ;;
    esac
done
outputs_dir="${eval_dir}/outputs/songs$suffix"
hashes_dir="${hashes_dir}/songs${suffix}"
mkdir -p "$hashes_dir"

if $generate; then
    for dir in "$outputs_dir"/*; do
        # calculate hash of playlist.m3u in each dir and save to hashes_dir
        if [ -d "$dir" ]; then
            bench=$(basename "$dir")
            md5sum "$dir/playlist.m3u" > "$hashes_dir/$bench.md5sum"
        fi
    done
    exit 0
fi

cd "$outputs_dir" || exit 1
status=0
for dir in *; do
    # check hash of playlist.m3u in each dir
    if [ -d "$dir" ]; then
        bench=$(basename "$dir")
        md5sum --check --quiet --status "$hashes_dir/$bench.md5sum"
        status=$?
    fi
done
echo playlist-creation $status