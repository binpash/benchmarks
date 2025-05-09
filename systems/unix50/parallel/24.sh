#!/bin/bash

# 9.2: extract the word BELL
# cat $1 | cut -c 1-1 | tr -d '\n'

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE" -k "cut -c 1-1" < "$1" | tr -d '\n'