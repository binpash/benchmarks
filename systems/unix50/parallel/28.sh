#!/bin/bash

# 9.6: Follow the directions for grep
# cat $1 | tr ' ' '\n' | grep '[A-Z]' | sed 1d | sed 3d | sed 3d | tr '[a-z]' '\n' | grep '[A-Z]' | sed 3d | tr -c '[A-Z]' '\n' | tr -d '\n'

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE" -k "tr ' ' '\n' | grep '[A-Z]'" < "$1" | sed 1d | sed 3d | sed 3d | tr '[a-z]' '\n' | grep '[A-Z]' | sed 3d | tr -c '[A-Z]' '\n' | tr -d '\n'
