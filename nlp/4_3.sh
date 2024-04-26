#!/bin/bash 
# tag: bigrams.sh
# set -e

# Bigrams (contrary to our version, this uses intermediary files)
IN=${IN:-$PWD/pg}
OUT=${OUT:-$PWD/output/4_3/}
ENTRIES=${ENTRIES:-10}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES})
do
    cat $IN/$input | tr -sc '[A-Z][a-z]' '[\012*]' > ${OUT}/${input}.input.words
    tail +2  ${OUT}/${input}.input.words > ${OUT}/${input}.input.nextwords
    paste ${OUT}/${input}.input.words ${OUT}/${input}.input.nextwords | sort | uniq -c > ${OUT}/${input}.input.bigrams
done

echo 'done';
