#!/bin/bash


IN=${IN:-$SUITE_DIR/inputs/pg}
INPUT2=${INPUT2:-$SUITE_DIR/inputs/exodus}
OUT=${1:-$SUITE_DIR/outputs/8.3_3/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

pure_func() {
    input=$1
    input2=$2
    paste <(tr -c 'A-Za-z' '[\n*]' < "$IN/$input" | sort -u) \
          <(tr -sc '[A-Z][a-z]' '[\012*]' < "$input2" | sort -u) |
        uniq -c | head
}

export -f pure_func

for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1); do
    pure_func "$input" "$INPUT2" > "${OUT}/${input}.out" &
done

wait

echo 'done'
