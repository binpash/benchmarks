#!/bin/bash

# 3.1: get lowercase first letter of last names (awk)
# cat $1 | cut -d ' ' -f 2 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE"  -k "cut -d ' ' -f 2 | cut -c 1-1" < "$1" | tr -d '\n' | tr '[A-Z]' '[a-z]'
