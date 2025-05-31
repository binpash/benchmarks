#!/bin/bash

eval_dir=$(realpath "$(dirname "$0")")
cd "$script_dir" || exit 1

export LC_ALL=C
suffix=""
generate=false
for arg in "$@"; do
    case "$arg" in
        --generate) generate=true ;;
        --small) suffix=".small" ;;
        --min) suffix=".min" ;;
    esac
done


hashes_dir="$eval_dir/hashes"

if $generate; then
    mkdir -p "$hashes_dir"

    echo "Generating file list for playlist-creation"

    outputs_dir="outputs/songs$suffix"
    
    filelist="$hashes_dir/songs$suffix.files"
    cd "$outputs_dir" || exit 1
    > "$filelist"
    for dir in *; do
        if [ -d "$dir" ] && [ -f "$dir/playlist.m3u" ]; then
            echo "$dir/playlist.m3u" >> "$filelist"
        fi
    done

    cd "$eval_dir" || exit 1
    echo "File list generated at $filelist"

    echo "Generating hashes for image-annotation"

    hashes_dir="${hashes_dir}/jpg$suffix"
    mkdir -p "$hashes_dir"

    outputs_dir="outputs/jpg$suffix"
    bench=image-annotation$suffix
    md5sum $outputs_dir/* > "$hashes_dir/$bench.md5sum"
    echo "Hashes generated at $hashes_dir/$bench.md5sum"
    exit 0
fi

outputs_dir="outputs/songs$suffix"
filelist="$hashes_dir/songs$suffix.files"

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

cd "$eval_dir" || exit 1

bench=image-annotation$suffix

hashes_dir="${hashes_dir}/jpg$suffix"
outputs_dir="outputs/jpg$suffix"

if [ ! -d "$outputs_dir" ]; then
    echo "Outputs directory not found: $outputs_dir"
    echo $bench 1
    exit 1
fi
if [ ! -d "$hashes_dir" ]; then
    echo "Hashes directory not found: $hashes_dir"
    echo $bench 1
    exit 1
fi
md5sum --check --quiet --status $hashes_dir/$bench.md5sum
echo $bench $?