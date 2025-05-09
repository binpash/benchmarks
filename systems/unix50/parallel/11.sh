#!/bin/bash

# 4.5: 4.4 + pawns
# cat $1 | tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P' | sort | uniq -c | sort -nr

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE"  -k "tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P'" < "$1" | sort | uniq -c | sort -nr