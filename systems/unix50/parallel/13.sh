#!/bin/bash

# 5.1: extract hello world
# cat $1 | grep 'print' | cut -d "\"" -f 2 | cut -c 1-12

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE" -k "grep 'print' | cut -d '\"' -f 2 | cut -c 1-12" < "$1"