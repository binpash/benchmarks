#!/bin/bash
# tag: count_vowel_seq
# set -e 

IN=${IN:-$PWD/pg/}
OUT=${OUT:-$PWD/output/2_2/}
ENTRIES=${ENTRIES:-10}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr 'a-z' '[A-Z]' | tr -sc 'AEIOU' '[\012*]'| sort | uniq -c  > ${OUT}/${input}.out
done

echo 'done';
