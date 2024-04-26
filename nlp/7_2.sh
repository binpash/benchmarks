#!/bin/bash
# set -e
# tag: count_consonant_sequences

IN=${IN:-$PWD/pg}
OUT=${OUT:-$PWD/output/7_2/}
ENTRIES=${ENTRIES:-10}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr '[a-z]' '[A-Z]' | tr -sc 'BCDFGHJKLMNPQRSTVWXYZ' '[\012*]' | sort | uniq -c > ${OUT}/${input}.out
done

echo 'done';
