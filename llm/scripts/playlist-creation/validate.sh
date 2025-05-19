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
filelist="$hashes_dir/songs$suffix.files"

if $generate; then
    mkdir -p "$hashes_dir"
    cd "$outputs_dir" || exit 1
    > "$filelist"
    for dir in *; do
        if [ -d "$dir" ] && [ -f "$dir/playlist.m3u" ]; then
            echo "$dir/playlist.m3u" >> "$filelist"
        fi
    done
    echo "Generated file list: $filelist"
    exit 0
fi

cd "$outputs_dir" || exit 1
status=0

if [ -f "$filelist" ]; then
    while read -r file; do
        if [ ! -f "$file" ]; then
            echo "File $file not found"
            status=1
        fi
    done < "$filelist"
else
    echo "File list not found: $filelist"
    status=1
fi

echo "playlist-creation $status"
exit $status
