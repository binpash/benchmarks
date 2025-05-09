#!/bin/sh 
# tag: sort_words_by_num_of_syllables

mkdir -p "$OUTPUT_DIR"

pure_func() {
    input=$1
    TEMPDIR=$(mktemp -d)
    cat > ${TEMPDIR}/${input}.words
    tr -sc '[AEIOUaeiou\012]' ' ' < ${TEMPDIR}/${input}.words | awk '{print NF}' > ${TEMPDIR}/${input}.syl
    paste ${TEMPDIR}/${input}.syl ${TEMPDIR}/${input}.words | sort -nr | sed 5q
    rm -rf ${TEMPDIR}
}

for input in $(ls $INPUT_FILE | xargs -I arg1 basename arg1)
do
    cat $INPUT_FILE/$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort -u | pure_func $input > $OUTPUT_DIR/$input.out
done
