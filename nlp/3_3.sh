#!/bin/bash
# tag: sort_words_by_rhyming.sh
# set -e

IN=${IN:-$PWD/pg}
OUT=${OUT:-$PWD/output/3_3/}
ENTRIES=${ENTRIES:-10}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | sort | uniq -c | rev | sort | rev > ${OUT}/${input}.out
done

echo 'done';
