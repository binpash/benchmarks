#!/bin/bash

# 4.6: piece used the most by Belle
# cat $1 | tr ' ' '\n' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P' | sort -r | uniq | head -n 3 | tail -n 1

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE" -k "tr ' ' '\n' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P'" < "$1" | sort -r | uniq | head -n 3 | tail -n 1

