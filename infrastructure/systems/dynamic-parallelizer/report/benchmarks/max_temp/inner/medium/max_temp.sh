#!/bin/bash

## Processing files and data per year
for year in $(seq $FROM $TO); do
    echo "Processing year: $year"
    find "$RESOURCE_DIR/$year" -type f -name '*.gz' |
    xargs -I {} gunzip -c {} > $RESOURCE_DIR/$year.txt

    ## Processing
    cat "$RESOURCE_DIR/$year.txt" |
    cut -c 89-92 |
    grep -v 999 |
    sort -rn |
    head -n1

    cat "$RESOURCE_DIR/$year.txt" |
    cut -c 89-92 |
    grep -v 999 |
    sort -n |
    head -n1

    cat "$RESOURCE_DIR/$year.txt" |
    cut -c 89-92 |
    grep -v 999 |
    awk "{ total += \$1; count++ } END { print total/count }"
done
