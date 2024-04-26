#!/bin/bash
# tag: resize image
IN=${IN:-inputs/jpg}
OUT=${OUT:-outputs/jpg}
mkdir -p ${OUT}
for i in $IN/*.jpg;
do
    out=$OUT/$(basename -- $i)
    convert -resize 70% "$i" "$out";
done

echo 'done';
