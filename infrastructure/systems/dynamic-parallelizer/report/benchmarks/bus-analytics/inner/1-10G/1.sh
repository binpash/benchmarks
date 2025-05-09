#!/bin/bash
# Vehicles on the road per day

cat $INPUT_FILE | sed 's/T..:..:..//' |
cut -d ',' -f 1,3 |
sort -u |
cut -d ',' -f 1 |
sort |
uniq -c |
awk -v OFS="\t" "{print \$2,\$1}"
