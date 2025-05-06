#!/bin/bash
# tag: wav-to-mp3
# inputs: $1=absolute source directory path with .wav's, $2=destination directory for output wavs

# Overwrite HOME variable
export HOME="$1"

mkdir -p $2

for i in ~/*;
do
    out="$2/$(basename $i).mp3"
    cat $i | ffmpeg -y -i pipe:0 -f mp3 -ab 192000 pipe:1 2>/dev/null > $out
done
