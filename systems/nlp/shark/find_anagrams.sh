#!/bin/bash 
# tag: find_anagrams.sh
# set -e

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/8.3_2/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

pure_func() {
    input=$1
    tr -c 'A-Za-z' '[\n*]' < "$IN/$input" | grep -v "^\s*$" | sort -u | rev | sort | uniq -c | awk '$1 >= 2 {print $2}'
}

export -f pure_func

for input in $(ls ${IN} | head -n ${ENTRIES} | xargs -I arg1 basename arg1)
do
    pure_func "$input" > "${OUT}/${input}.out" &
done

wait

echo 'done';
# rm -rf "$OUT"
