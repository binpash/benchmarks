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
mkdir -p "$hashes_dir"

hashfile="$hashes_dir/songs$suffix.md5sum"
export LC_ALL=C

if $generate; then
    mkdir -p "$hashes_dir"
    cd "$outputs_dir" || exit 1
    > "$hashfile"
    for dir in *; do
        if [ -d "$dir" ]; then
            playlist_path="$dir/playlist.m3u"
            if [ -f "$playlist_path" ]; then
                md5sum "$playlist_path" >> "$hashfile"
            fi
        fi
    done
    exit 0
fi

cd "$outputs_dir" || exit 1
status=0

if [ -f "$hashfile" ]; then
    md5sum -c "$hashfile" --quiet || status=1
else
    echo "Hash file not found: $hashfile"
    status=1
fi

echo "playlist-creation $status"