#!/bin/bash

# 7.1: identify number of AT&T unix versions
# cat $1 | cut -f 1 | grep 'AT&T' | wc -l

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE" -k "cut -f 1 | grep 'AT&T'" < "$1" | wc -l