#!/bin/sh
# tag: compare_exodus_genesis.sh

INPUT2="$INPUT_TOP/exodus"

pure_func() {
    input=$1
    input2=$2
    TEMPDIR=$(mktemp -d)
    cat > ${TEMPDIR}/${input}1.types
    cat  ${input2} | tr -sc '[A-Z][a-z]' '[\012*]' | sort -u > ${TEMPDIR}/${input}2.types
    sort ${TEMPDIR}/${input}1.types ${TEMPDIR}/${input}2.types ${TEMPDIR}/${input}2.types | uniq -c | head 
    rm -rf ${TEMPDIR}
}

for input in $(ls $INPUT_FILE | xargs -I arg1 basename arg1)
do
    cat $INPUT_FILE/$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort -u | pure_func $input $INPUT2 > $OUTPUT_DIR/$input.out
done
