#!/bin/bash
# tag: merge_upper
# set -e

# Merge upper and lower counts
IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/2_1/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

# for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1)
# do
#     cat $IN/$input | tr '[a-z]' '[A-Z]' |  tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort | uniq -c > ${OUT}/${input}.out
# done

for input in $(ls "$IN" | head -n "$ENTRIES"); do
    tr '[a-z]' '[A-Z]' < "$IN/$input" | \
    tr -c 'A-Za-z' '[\n*]' | \
    grep -v "^\s*$" | \
    sort | uniq -c > "${OUT}/${input}.out" &
done
wait

echo 'done';
# rm -rf "$OUT"
