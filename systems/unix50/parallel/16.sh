#!/bin/bash

# 7.2: find  most frequently occurring machine
# cat $1 | cut -f 2 | sort -n | uniq -c | sort -nr | head -n 1 | tr -s ' ' '\n' | tail -n 1

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE" -k "cut -f 2" < "$1" | sort | uniq -c | sort -nr | head -n 1 | tr -s ' ' '\n' | tail -n 1