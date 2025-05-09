#!/bin/bash
# tag: count_words

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/1_1/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1)
do
    tr -c 'A-Za-z' '[\n*]' < "$IN/$input" | grep -v "^\s*$" | sort | uniq -c > "${OUT}/${input}.out" &
done

wait

echo 'done';