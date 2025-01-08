#!/bin/bash

# 11.1: year Ritchie and Thompson receive the Hamming medal
# cat $1 | grep 'UNIX' | cut -f 1

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE" -k "grep 'UNIX' | cut -f 1" < "$1"