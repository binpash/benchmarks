#!/bin/bash
# Hours monitored each day

# diff out{1,}

sed 's/T\(..\):..:../,\1/' "$1" |  # keep times only
  cut -d ',' -f 1,2 |           # keep only time and date
  sort -u |                     # removing duplicate entries
  cut -d ',' -f 1 |             # keep only date
  sort |                        # preparing for uniq
  uniq -c |                     # count unique dates
  awk '{print $2, $1}'         # print first date, then count