#!/bin/bash

# 8.3: find names of the four people most involved with unix
# cat $1 | grep '(' | cut -d '(' -f 2 | cut -d ')' -f 1 | head -n 1

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE" -k "grep '(' | cut -d '(' -f 2 | cut -d ')' -f 1" < "$1" | head -n 1