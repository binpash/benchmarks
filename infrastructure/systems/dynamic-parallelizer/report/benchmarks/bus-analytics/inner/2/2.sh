#!/bin/bash
# Days a vehicle is on the road

cat $INPUT_FILE |
  sed 's/T..:..:..//' |
  cut -d ',' -f 3,1 |
  sort -u |
  cut -d ',' -f 2 |
  sort |
  uniq -c |
  sort -k1n |
  awk -v OFS="\t" "{print \$2,\$1}"
