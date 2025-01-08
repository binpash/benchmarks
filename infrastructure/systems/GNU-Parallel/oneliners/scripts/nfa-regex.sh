#!/bin/bash
# Match complex regular-expression over input

# cat $1 | tr A-Z a-z | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4'

# using GNU parallel

cat "$1" | parallel --pipe -k --block 100K "tr A-Z a-z | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4'"
