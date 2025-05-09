#!/bin/sh 
# tag: find_anagrams.sh

mkdir -p "$OUTPUT_DIR"

pure_func() {
    input=$1
    TEMPDIR=$(mktemp -d)
    sort -u > ${TEMPDIR}/${input}.types
    rev < ${TEMPDIR}/${input}.types > ${TEMPDIR}/${input}.types.rev
    sort ${TEMPDIR}/${input}.types ${TEMPDIR}/${input}.types.rev | uniq -c | awk "\$1 >= 2 {print \$2}"
    rm -rf ${TEMPDIR}
}

for input in $(ls $INPUT_FILE | xargs -I arg1 basename arg1)
do
    cat $INPUT_FILE/$input |  tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | pure_func $input > $OUTPUT_DIR/$input.out
done
