#!/bin/bash
# Calculate sort twice

IN=${IN:-./input_txt/100M.txt}

cat $IN | tr A-Z a-z | sort | sort -r
