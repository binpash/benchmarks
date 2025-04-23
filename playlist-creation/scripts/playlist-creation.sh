#!/bin/bash
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/playlist-creation"
inputs_dir="$eval_dir/min_inputs"
outputs_dir="$eval_dir/outputs"

IN=${1:-"$inputs_dir"}
OUT="${2:-"$outputs_dir"}"
mkdir -p "$OUT"

files=$(find "$IN" -type f -name "*.wav" | sort)
num_files=$(printf '%s\n' "$files" | wc -l)

llm embed-multi -m clap songs --binary --files "$IN" '*.wav'

first_song=$(printf '%s\n' "$files" | head -n 1)
last_song=$(printf '%s\n' "$files" | tail -n 1)

playlist_path="$OUT/playlist.m3u"
llm interpolate songs "$first_song" "$last_song" -n "$num_files" | jq .[] > "$outputs_dir/$playlist_path"
