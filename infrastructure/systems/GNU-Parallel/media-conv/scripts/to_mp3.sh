#!/bin/bash
# tag: wav-to-mp3
# inputs: $1=absolute source directory path with .wav's, $2=destination directory for output wavs

# mkdir -p $2

# pure_func(){
#     ffmpeg -y -i pipe:0 -f mp3 -ab 192000 pipe:1  2>/dev/null
# }
# export -f pure_func

# for i in $1/*;
# do
#     out="$2/$(basename $i).mp3"
#     cat $i | pure_func > $out
# done

# Using GNU parallel:

mkdir -p "$2"

pure_func(){
    ffmpeg -y -i pipe:0 -f mp3 -ab 192000 pipe:1 2>/dev/null
}
export -f pure_func

export DEST_DIR="$2"

parallel "cat {} | pure_func > $DEST_DIR/{/}.mp3" ::: "$1"/*