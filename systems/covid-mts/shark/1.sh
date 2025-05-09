#!/bin/bash
# Vehicles on the road per day

sed 's/T..:..:..//' "$1" |        # hide times directly from the input file
  cut -d ',' -f 1,3 |             # keep only day and bus no
  sort -u |                       # remove duplicate records
  cut -d ',' -f 1 |               # keep all dates
  sort |                          # prepare for uniq
  uniq -c |                       # count unique dates
  awk '{print $2, $1}'            # print date first, then count

# diff out{1,}
