#!/bin/bash

mkdir -p "$2"

pure_func() {
    ffmpeg -y -i pipe:0 -f mp3 -ab 192000 pipe:1 2>/dev/null
}
export -f pure_func

export dest="$2"

find "$1" -type f | parallel ' \
    out="$dest/$(basename {}).mp3"; \
    cat {} | pure_func > "$out" \
'