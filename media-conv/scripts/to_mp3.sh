#!/bin/bash
# tag: wav-to-mp3
# inputs: $1=absolute source directory path with .wav's, $2=destination directory for output wavs

# Overwrite HOME variable
export HOME="$1"

mkdir -p $2

for i in ~/*;
do
    out="$2/$(basename "$i").mp3"
    cat "$i" |
        ffmpeg -hide_banner -loglevel error \
        -y -i pipe:0 \
        -acodec libmp3lame -b:a 192k \
        -ar 44100 -ac 2 \
        -fflags +bitexact \
        -flags:a +bitexact \
        -write_xing 0 \
        -id3v2_version 0 \
        -map_metadata -1 \
        -f mp3 pipe:1 \
        > "$out" 2>/dev/null
done
