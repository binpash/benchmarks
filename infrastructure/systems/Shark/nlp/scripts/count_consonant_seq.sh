#!/bin/bash
# set -e
# tag: count_consonant_sequences

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/7_2/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

# for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1)
# do
#     cat $IN/$input | tr '[a-z]' '[A-Z]' | tr -sc 'BCDFGHJKLMNPQRSTVWXYZ' '[\012*]' | sort | uniq -c > ${OUT}/${input}.out
# done

for input in $(ls "$IN" | head -n "$ENTRIES"); do
    tr '[a-z]' '[A-Z]' < "$IN/$input" | \
    tr -sc 'BCDFGHJKLMNPQRSTVWXYZ' '[\012*]' | \
    sort | uniq -c > "${OUT}/${input}.out" &
done

wait

echo 'done';
# rm -rf ${OUT}
