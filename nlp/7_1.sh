#!/bin/bash
# tag: count_morphs
# set -e

IN=${IN:-$PWD/pg}
OUT=${OUT:-$PWD/output/7_1/}
ENTRIES=${ENTRIES:-10}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | sed 's/ly$/-ly/g' | sed 's/ .*//g' | sort | uniq -c > ${OUT}/${input}.out
done

echo 'done';
