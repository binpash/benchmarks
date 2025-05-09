#!/bin/bash
# compress all files in a directory

# Using GNU parallel:

mkdir -p "$2"

compress_file() {
  input_file="$1"
  output_file="$2/$(basename "$input_file").zip"
  gzip --no-name -c "$input_file" > "$output_file"
}

export -f compress_file

find "$1" -type f -name "*.pcapng" | parallel --jobs "$(nproc)" compress_file {} "$2"