#!/bin/bash
# tag: count_words

IN=${IN:-$PWD/pg/}
OUT=${OUT:-$PWD/output/1_1/}
ENTRIES=${ENTRIES:-10}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | sort | uniq -c > ${OUT}/${input}.out
done

echo 'done';
