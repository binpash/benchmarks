#!/bin/bash

# 9.3: animal that used to decorate the Unix room
# cat $1 | cut -c 1-2 | tr -d '\n'

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE" -k "cut -c 1-2" < "$1" | tr -d '\n'