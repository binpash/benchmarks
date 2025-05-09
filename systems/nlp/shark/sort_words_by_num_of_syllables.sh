#!/bin/bash 
# tag: sort_words_by_num_of_syllables
# set -e

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/8.1/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"


pure_func() {
    input=$1
    tr -sc '[AEIOUaeiou\012]' ' ' < "$IN/$input" | awk '{print NF}' |
    paste - <(tr -c 'A-Za-z' '[\n*]' < "$IN/$input" | sort -u) | sort -nr | sed 5q
}
for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1); do
    pure_func "$input" > "${OUT}/${input}.out" &
done
wait
echo 'done';
# rm -rf "$OUT"
