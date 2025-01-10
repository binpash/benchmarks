#!/bin/bash
# Calculate sort twice

# cat $1 | tr A-Z a-z | sort | sort -r

cat "$1" | parallel -k --pipe --block "$BLOCK_SIZE" "tr A-Z a-z" | sort | sort -r
