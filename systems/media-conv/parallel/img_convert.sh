#!/bin/bash

# Using GNU parallel:

mkdir -p "$2"

pure_func () {
    convert -resize 70% "-" "-"
}
export -f pure_func

export dest_dir="$2"

find "$1" -type f | parallel --jobs "$(nproc)" \
    'cat {} | pure_func > "$dest_dir/{/}"'