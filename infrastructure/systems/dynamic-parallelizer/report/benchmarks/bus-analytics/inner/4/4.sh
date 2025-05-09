#!/bin/bash
# Hours monitored each day

cat $INPUT_FILE |
  sed 's/T\(..\):..:../,\1/' |
  cut -d ',' -f 1,2 |
  sort -u |
  cut -d ',' -f 1 |
  sort |
  uniq -c |
  awk -v OFS="\t" "{print \$2,\$1}"
