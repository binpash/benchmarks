#!/bin/bash 

IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/4_3/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

pure_func() {
    input=$1
    IN_DIR=$2
    OUT_DIR=$3
    TEMPDIR=$(mktemp -d)

    # Process the input file
    cat "$IN_DIR/$input" | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" > ${TEMPDIR}/${input}.input.words
    tail +2 ${TEMPDIR}/${input}.input.words > ${TEMPDIR}/${input}.input.nextwords
    paste ${TEMPDIR}/${input}.input.words ${TEMPDIR}/${input}.input.nextwords | sort | uniq -c > "$OUT_DIR/${input}.input.bigrams.out"

    # Cleanup temporary directory
    rm -rf ${TEMPDIR}
}
export -f pure_func

ls "${IN}" | head -n "${ENTRIES}" | parallel -j "$(nproc)" pure_func {} "${IN}" "${OUT}"

echo "done"