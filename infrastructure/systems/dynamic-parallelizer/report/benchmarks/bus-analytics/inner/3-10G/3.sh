#!/bin/bash
# Hours each vehicle is on the road

cat $INPUT_FILE |
  sed 's/T\(..\):..:../,\1/' |
  cut -d ',' -f 1,2,4 |
  sort -u |
  cut -d ',' -f 3 |
  sort |
  uniq -c |
  sort -k1n |
  awk -v OFS="\t" "{print \$2,\$1}"