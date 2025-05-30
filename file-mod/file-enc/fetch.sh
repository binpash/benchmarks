#!/bin/bash

TOP=$(git rev-parse --show-toplevel)

input_dir="${TOP}/file-enc/inputs"
mkdir -p "$input_dir"

URL="https://atlas.cs.brown.edu/data"
ZIP_DST="$input_dir/pcaps.zip"

size=full
for arg in "$@"; do
    case "$arg" in
        --small) size=small ;;
        --min)   size=min ;;
    esac
done

if [ -d "$input_dir/pcaps_$size" ]; then
    echo "Directory $input_dir/pcaps_$size already exists. Skipping download."
    exit 0
fi

if [ "$size" = "min" ]; then
    mkdir -p "$input_dir/pcaps_$size"
    cp min_inputs/* "$input_dir/pcaps_$size/"
    exit 0
fi

if [ "$size" = "small" ]; then
    wget --no-check-certificate $URL/pcaps.zip -O "$input_dir/pcaps_$size.zip"
    unzip "$input_dir/pcaps_$size.zip" -d "$input_dir"
    mv "$input_dir/pcaps/" "$input_dir/pcaps_$size/"
    rm "$input_dir/pcaps_$size.zip"
    exit 0
fi

wget --no-check-certificate $URL/pcaps.zip -O "$input_dir/pcaps.zip"
unzip "$input_dir/pcaps.zip" -d "$input_dir"
mv "$input_dir/pcaps/" "$input_dir/pcaps_$size/"
rm "$input_dir/pcaps.zip"

wget --no-check-certificate "$URL"/pcaps_large.zip -O "$ZIP_DST"
unzip "$ZIP_DST" -d "$input_dir/pcaps_$size"
rm "$ZIP_DST"
