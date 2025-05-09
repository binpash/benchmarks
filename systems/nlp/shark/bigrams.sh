#!/bin/bash 

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/4_3/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

pure_func() {
    input=$1

    # Generate bigrams directly in memory
    tr -c 'A-Za-z' '[\n*]' < "$IN/$input" | grep -v "^\s*$" | \
    paste - <(tr -c 'A-Za-z' '[\n*]' < "$IN/$input" | grep -v "^\s*$" | tail -n +2)
}

export -f pure_func

# Process each input file
for input in $(find "$IN" -type f | head -n ${ENTRIES}); do
    pure_func "$(basename "$input")" | \
    sort | uniq -c > "${OUT}/$(basename "$input").input.bigrams.out" &
done

# Wait for all background processes to finish
wait

echo 'done'
