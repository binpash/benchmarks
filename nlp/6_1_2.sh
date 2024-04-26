#!/bin/bash
# tag: uppercase_by_type
# set -e

IN=${IN:-$PWD/pg}
OUT=${OUT:-$PWD/output/6_1_2/}
ENTRIES=${ENTRIES:-10}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u | grep -c '^[A-Z]' > ${OUT}/${input}.out
done

echo 'done';
