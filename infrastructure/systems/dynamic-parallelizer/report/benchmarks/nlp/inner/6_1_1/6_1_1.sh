#!/bin/sh
# tag: uppercase_by_token

mkdir -p "$OUTPUT_DIR"

for input in $(ls $INPUT_FILE | xargs -I arg1 basename arg1)
do
    cat $INPUT_FILE/$input |  tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" | grep -c '^[A-Z]' > $OUTPUT_DIR/$input.out
done
