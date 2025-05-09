#!/bin/bash
# tag: bigrams_appear_twice.sh

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/8.2_2/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

pure_func() {
    local input="$1"
    tr -c 'A-Za-z' '[\n*]' < "$IN/$input" | grep -v "^\s*$" | \
    awk '{if (NR > 1) print prev, $0; prev = $0}' | \
    sort | uniq -c | \
    awk '$1 == 2 {print $2, $3}'
}

export -f pure_func

for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1)
do
    pure_func "$input" > "${OUT}/${input}.out" &
done

wait
echo 'done'
