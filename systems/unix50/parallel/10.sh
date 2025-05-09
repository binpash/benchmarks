#!/bin/bash

# 4.4: histogram of Belle's captures (-pawns) by each type of piece
# cat $1 | tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | grep '[KQRBN]' | cut -c 1-1 | sort | uniq -c | sort -nr

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE"  -k "tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | grep '[KQRBN]' | cut -c 1-1" < "$1" | sort | uniq -c | sort -nr