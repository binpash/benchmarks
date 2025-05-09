#!/bin/bash

# 8.4: find longest words without hyphens
# cat $1 | tr -c "[a-z][A-Z]" '\n' | sort | awk "length >= 16"

# Using GNU parallel:
parallel --jobs "$jobs" --pipe --block "$BLOCK_SIZE" -k "tr -c '[a-z][A-Z]' '\n'" < "$1" | sort | awk 'length >= 16'
