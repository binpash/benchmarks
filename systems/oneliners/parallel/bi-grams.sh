#!/bin/bash
# Find all 2-grams in a piece of text

. ./scripts/bi-gram.aux.sh
cat $1 |
  tr -c 'A-Za-z' '[\n*]' | 
  grep -v "^\s*$" |
  tr A-Z a-z |
  bigrams_aux |
  sort |
  uniq

# No opportunity for GNU parallel