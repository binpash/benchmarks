#!/bin/bash

# 9.4: four corners with E centered, for an "X" configuration
# cat $1 | tr ' ' '\n' | grep "\"" | sed 4d | cut -d "\"" -f 2 | tr -d '\n'

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE" -k "tr ' ' '\n' | grep '\"' | cut -d '\"' -f 2" < "$1" | sed 4d | tr -d '\n'
