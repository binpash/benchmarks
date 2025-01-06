#!/bin/bash
# Hours monitored each day

# <in.csv sed 's/T\(..\):..:../,\1/' |
# awk -F, '!seen[$1 $2] {hours[$1]++; seen[$1 $2] = 1}
#    END { OFS = "\t"; for (d in hours) print d, hours[d]}' |
#   sort

# curl https://balab.aueb.gr/~dds/oasa-$(date --date='1 days ago' +'%y-%m-%d').bz2 |
#   bzip2 -d |                  # decompress
# Replace the line below with the two lines above to stream the latest file
# cat "$1" |                        # assumes saved input
#   sed 's/T\(..\):..:../,\1/' |  # keep times only
#   cut -d ',' -f 1,2 |           # keep only time and date
#   sort -u |                     # removing duplicate entries
#   cut -d ',' -f 1 |             # keep only date
#   sort |                        # preparing for uniq
#   uniq -c |                     # count unique dates
#   awk "{print \$2,\$1}"         # print first date, then count

# diff out{1,}

# Using GNU parallel:

INPUT="$1"

process_chunk() {
  local chunk="$1"
  sed 's/T\(..\):..:../,\1/' "$chunk" |  # keep times only
  cut -d ',' -f 1,2                      # keep only time and date
}
export -f process_chunk

lines=$(wc -l < "$1")
nproc=$(nproc)
chunk_size=$((lines / nproc))


split -l "$chunk_size""$INPUT" chunk_

ls chunk_* | parallel -j "$(nproc)" process_chunk > combined.tmp

cat combined.tmp |
  sort -u |                       # global deduplication
  cut -d ',' -f 1 |               # keep only date
  sort |                          # prepare for counting
  uniq -c |                       # count unique dates
  awk '{print $2,$1}'             # print date, then count

rm chunk_*
rm combined.tmp