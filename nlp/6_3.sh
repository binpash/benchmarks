#!/bin/bash
# tag: words_no_vowels
# set -e

IN=${IN:-$PWD/pg}
OUT=${OUT:-$PWD/output/6_3/}
ENTRIES=${ENTRIES:-10}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | grep -vi '[aeiou]' | sort | uniq -c > ${OUT}/${input}.out
done

echo 'done';
