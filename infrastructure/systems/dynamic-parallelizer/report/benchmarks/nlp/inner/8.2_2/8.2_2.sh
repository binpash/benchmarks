#!/bin/sh 
# tag: bigrams_appear_twice.sh

mkdir -p "$OUTPUT_DIR"

pure_func() {
    input=$1
    TEMPDIR=$(mktemp -d)
    cat > ${TEMPDIR}/${input}.input.words
    tail +2 ${TEMPDIR}/${input}.input.words > ${TEMPDIR}/${input}.input.nextwords
    paste ${TEMPDIR}/${input}.input.words ${TEMPDIR}/${input}.input.nextwords | sort | uniq -c > ${TEMPDIR}/${input}.input.bigrams
    awk "\$1 == 2 {print \$2, \$3}" ${TEMPDIR}/${input}.input.bigrams
    rm -rf {TEMPDIR}
}

for input in $(ls $INPUT_FILE | xargs -I arg1 basename arg1)
do
    cat $INPUT_FILE/$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | pure_func $input > $OUTPUT_DIR/$input.out
done
