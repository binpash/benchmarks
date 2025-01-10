#!/bin/bash
# Match complex regular-expression over input

# cat $1 | tr A-Z a-z | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4'

tr A-Z a-z < $1 | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4'