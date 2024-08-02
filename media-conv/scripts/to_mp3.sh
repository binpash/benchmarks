#!/bin/bash
# tag: wav-to-mp3
# IN=${IN:-/media-conv/wav}
# OUT=${OUT:-$DISH_TOP/evaluation/media-conv/outputs/mp3}
mkdir -p $2

pure_func(){
    ffmpeg -y -i pipe:0 -f mp3 -ab 192000 pipe:1  2>/dev/null
}
export -f pure_func

for item in $(ls $1);
do
    pkg_count=$((pkg_count + 1));
    out="$2/$(basename $item).mp3"
    cat $item | pure_func > $out
done

echo 'done';
