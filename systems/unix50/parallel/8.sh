#!/bin/bash

# 4.2: find pieces captured by Belle
# cat $1 | tr ' ' '\n' | grep 'x' | grep '\.' | wc -l

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE"  -k "tr ' ' '\n' | grep 'x' | grep '\.'" < "$1" | wc -l