#!/bin/bash
# tag: uppercase_by_token
# set -e

IN=${IN:-$PWD/pg}
OUT=${OUT:-$PWD/output/6_1_1/}
ENTRIES=${ENTRIES:-10}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | grep -c '^[A-Z]' > ${OUT}/${input}.out
done

echo 'done';
