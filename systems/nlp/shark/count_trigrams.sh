#!/bin/bash 

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/4_3b/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

pure_func() {
    input=$1
    tr -c 'A-Za-z' '[\n*]' < "$IN/$input" | grep -v "^\s*$" |
    tail -n +2 | paste - <(tail -n +2 "$IN/$input" | tail -n +2) | sort | uniq -c
}

export -f pure_func

# Parallelize the processing
for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1); do
    pure_func "$input" > "${OUT}/${input}.trigrams" &
done

wait

echo 'done'

# rm -rf ${OUT}
