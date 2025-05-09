#!/bin/bash

# 7.3: all the decades in which a unix version was released
# cat $1 | cut -f 4 | sort -n | cut -c 3-3 | uniq | sed s/\$/'0s'/

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE" -k "cut -f 4" < "$1" | sort -n | cut -c 3-3 | uniq | sed s/\$/'0s'/