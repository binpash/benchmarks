#!/bin/bash
# Match complex regular-expression over input

IN=${IN:-./input_txt/10M.txt}

cat $IN | tr A-Z a-z | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4'
