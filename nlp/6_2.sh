#!/bin/bash
# tag: four-letter words
# set -e

# the original script has both versions
IN=${IN:-$PWD/pg}
OUT=${OUT:-$PWD/output/6_2/}
ENTRIES=${ENTRIES:-10}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | grep -c '^....$' > ${OUT}/${input}.out0
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u | grep -c '^....$'  > ${OUT}/${input}.out1
done

echo 'done';
