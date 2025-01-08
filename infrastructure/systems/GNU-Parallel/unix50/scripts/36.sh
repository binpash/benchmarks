#!/bin/bash

# 11.2: most repeated first name in the list?
# cat $1 | cut -f 2 | cut -d ' ' -f 1 | sort | uniq -c | sort -nr | head -n 1 | fmt -w1 | sed 1d

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE" -k "cut -f 2 | cut -d ' ' -f 1" < "$1" | sort | uniq -c | sort -nr | head -n 1 | fmt -w1 | sed 1d