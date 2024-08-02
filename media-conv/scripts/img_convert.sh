#!/bin/bash
# tag: resize image 
# IN=${JPG:-/media-conv/jpg}
# OUT=${OUT:-$DISH_TOP/evaluation/media-conv/outputs/jpg}
mkdir -p $2

pure_func () {
     convert -resize 70% "-" "-"
}
export -f pure_func

for i in $(ls $1/*.jpg); 
do 
    out=$2/$(basename -- $i)
    cat $i | pure_func > $out; 
done

echo 'done';
