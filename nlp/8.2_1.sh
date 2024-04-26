#!/bin/bash
# tag: vowel_sequences_gr_1K.sh
# set -e

IN=${IN:-$PWD/pg}
OUT=${OUT:-$PWD/output/8_2_1/}
ENTRIES=${ENTRIES:-10}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' | tr -sc 'AEIOUaeiou' '[\012*]' | sort | uniq -c | awk "\$1 >= 1000" > ${OUT}/${input}.out
done

echo 'done';
