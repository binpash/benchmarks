#!/bin/bash
# source: posh benchmark suite

# Overwrite HOME variable
export HOME="$1"
dest="$2"

mkdir -p "$dest"

resize_image () {
    mogrify -format gif -path "$dest" -thumbnail 100x100 "$1"
}

export -f resize_image

for file in ~/*; do
    resize_image "$file"
done
