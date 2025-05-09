#!/bin/sh
# tag: sort_words_by_rhyming.sh

mkdir -p "$OUTPUT_DIR"

for input in $(ls $INPUT_FILE | xargs -I arg1 basename arg1)
do
    cat $INPUT_FILE/$input | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | sort | uniq -c | rev | sort | rev > $OUTPUT_DIR/$input.out
done
