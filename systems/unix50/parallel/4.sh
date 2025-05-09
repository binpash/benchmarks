#!/bin/bash

# 1.3: sort top first names
# cat $1 | cut -d ' ' -f 1 | sort | uniq -c | sort -r

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE"  -k "cut -d ' ' -f 1" < "$1" | sort | uniq -c | sort -r