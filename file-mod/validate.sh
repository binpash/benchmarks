#!/bin/bash

cd "$(realpath "$(dirname "$0")")" || exit 1

outputs_dir="outputs"
hashes_dir="hashes"
TOP="$(git rev-parse --show-toplevel)"
eval_dir="${TOP}/file-mod"
suffix=".full"

mkdir -p "$hashes_dir"

generate=false
for arg in "$@"; do
    if [[ "$arg" == "--generate" ]]; then
        generate=true
        continue
    fi
    case "$arg" in
        --small) suffix=".small" ;;
        --min) suffix=".min" ;;
    esac
done

hash_audio_dir() {
    local src_dir=$1
    for src in $src_dir/*; do
        got_hash=$(ffmpeg -i "$src" -map 0:a -f md5 - 2>/dev/null)
        echo $got_hash $(realpath "--relative-to=$src_dir" "$src")
    done
}

if $generate; then
    md5sum $outputs_dir/compress_files$suffix/* > "$hashes_dir/compress_files$suffix.md5sum"
    md5sum $outputs_dir/encrypt_files$suffix/* > "$hashes_dir/encrypt_files$suffix.md5sum"
    md5sum $outputs_dir/img_convert$suffix/* > "$hashes_dir/img_convert$suffix.md5sum"
    md5sum $outputs_dir/thumbnail$suffix/* > "$hashes_dir/thumbnail$suffix.md5sum"
    hash_audio_dir "$eval_dir/$outputs_dir/to_mp3$suffix" > "$eval_dir/$hashes_dir/to_mp3$suffix.md5sum"
    echo "Generated hashes in $hashes_dir"
    exit 0
fi

status=0
if ! md5sum --check --quiet "$hashes_dir/encrypt_files$suffix.md5sum"; then
    status=1
fi
echo "encrypt_files $status"

status=0
if ! md5sum --check --quiet "$hashes_dir/compress_files$suffix.md5sum"; then
    status=1
fi
echo "compress_files $status"

status=0
if ! md5sum --check --quiet "$hashes_dir/img_convert$suffix.md5sum"; then
    status=1
fi
echo "img_convert $status"

status=0
if ! md5sum --check --quiet "$hashes_dir/thumbnail$suffix.md5sum"; then
    status=1
fi
echo "thumbnail $status"

hash_audio_dir "$eval_dir/$outputs_dir/to_mp3$suffix" | diff -q "$eval_dir/$hashes_dir/to_mp3$suffix.md5sum" -
echo "to_mp3 $?"
