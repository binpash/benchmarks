#!/bin/bash
# tag: trigram_rec
# set -e

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/6_1/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"


pure_func() {
    input=$1
    tr -sc '[AEIOUaeiou\012]' ' ' < "$IN/$input" | awk '{print NF}' |
    paste - <(tr -c 'A-Za-z' '[\n*]' < "$IN/$input" | sort -u) | sort -nr | sed 5q
}

export -f pure_func

for input in $(ls "$IN" | head -n "$ENTRIES"); do
    pure_func "$input" > "${OUT}/${input}.out" &
done

wait

echo 'done';
# rm -rf "$OUT"
