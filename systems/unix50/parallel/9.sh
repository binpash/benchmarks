#!/bin/bash

# 4.3: find pieces captured by Belle with a pawn
# cat $1 | tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | grep -v '[KQRBN]' | wc -l

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE"  -k "tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | grep -v '[KQRBN]'" < "$1" | wc -l