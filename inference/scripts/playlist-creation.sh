#!/bin/bash

IN="$1"
OUT="$2"
mkdir -p "$OUT"

# files=$(find "$IN" -type f \( -iname "*.mp3" -o -iname "*.wav" \) | sort)
# num_files=$(printf '%s\n' "$files" | wc -l)

for dir_path in "$IN"/*; do
    [ -d "$dir_path" ] || continue
    dir=$(basename "$dir_path")
    echo "Processing directory: $dir"

    files=$(find "$dir_path" -type f -name "*.mp3" | sort)
    num_files=$(printf '%s\n' "$files" | wc -l)

    abs_prefix=$(realpath "$dir_path")
    llm embed-multi -m clap songs --binary --files "$abs_prefix" '*.mp3' --prefix "$abs_prefix/"
    
    first_song=$(printf '%s\n' "$files" | head -n 1)
    last_song=$(printf '%s\n' "$files" | tail -n 1)

    mkdir -p "$OUT/$dir"
    playlist_path="$OUT/$dir/playlist.m3u"

    llm interpolate songs "$first_song" "$last_song" -n "$num_files" | jq .[] > "$playlist_path"
done