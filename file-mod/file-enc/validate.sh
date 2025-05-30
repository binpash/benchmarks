#!/bin/bash

# The output files .enc and .zip will vary across separate runs
# even for the same input file. Therefore, correctness check is omitted.

#md5sum compress_files.full/* > /benchmarks/file-enc/hashes/compress_files.full.md5sum
#md5sum encrypt_files.full/* > /benchmarks/file-enc/hashes/encrypt_files.full.md5sum
#md5sum compress_files.full/* > /benchmarks/file-enc/hashes/compress_files.small.md5sum
#md5sum encrypt_files.full/* > /benchmarks/file-enc/hashes/encrypt_files.small.md5sum

cd "$(realpath "$(dirname "$0")")" || exit 1

outputs_dir="outputs"
hashes_dir="hashes"

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

if $generate; then
    md5sum $outputs_dir/compress_files$suffix/* > "$hashes_dir/compress_files$suffix.md5sum"
    md5sum $outputs_dir/encrypt_files$suffix/* > "$hashes_dir/encrypt_files$suffix.md5sum"
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
