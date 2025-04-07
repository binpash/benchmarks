#!/bin/bash

src="$1"
dest="$2"

mkdir -p "$dest"

resize_image () {
    mogrify -format gif -path "$dest" -thumbnail 100x100 "$1"
}

export -f resize_image

for file in "$src"/*; do
    resize_image "$file"
done
