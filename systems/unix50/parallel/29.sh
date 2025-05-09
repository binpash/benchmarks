#!/bin/bash

# 9.7: Four corners
# cat $1 | sed 2d | sed 2d | tr -c '[A-Z]' '\n' | tr -d '\n'

# Using GNU parallel:
cat $1 | sed 2d | sed 2d | parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE" -k "tr -c '[A-Z]' '\n'" | tr -d '\n'
