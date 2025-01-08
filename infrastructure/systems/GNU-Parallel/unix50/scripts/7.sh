#!/bin/bash

# 4.1: find number of rounds
# cat $1 | tr ' ' '\n' | grep '\.' | wc -l

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE"  -k "tr ' ' '\n' | grep '\.'" < "$1" | wc -l