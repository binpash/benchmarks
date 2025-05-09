#!/bin/bash

# 10.3: extract Ritchie's username
# cat $1 | grep 'Bell' | cut -f 2 | head -n 1 | fmt -w1 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE" -k "grep 'Bell' | cut -f 2" < "$1" | head -n 1 | fmt -w1 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'