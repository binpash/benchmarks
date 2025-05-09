#!/bin/sh
# tag: count_morphs

mkdir -p "$OUTPUT_DIR"

for input in $(ls $INPUT_FILE | xargs -I arg1 basename arg1)
do
    cat $INPUT_FILE/$input | sed 's/ly$/-ly/g' | sed 's/ .*//g' | sort | uniq -c > $OUTPUT_DIR/$input.out
done
