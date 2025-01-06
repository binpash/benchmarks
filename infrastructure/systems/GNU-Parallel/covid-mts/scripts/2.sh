#!/bin/bash
# Days a vehicle is on the road

# <in.csv sed 's/T..:..:..//' |
# awk -F, '!seen[$1 $3] {onroad[$3]++; seen[$1 $3] = 1}
#    END { OFS = "\t"; for (d in onroad) print d, onroad[d]}' |
# sort -k2n >out1

# curl https://balab.aueb.gr/~dds/oasa-$(date --date='1 days ago' +'%y-%m-%d').bz2 |
#   bzip2 -d |                  # decompress
# Replace the line below with the two lines above to stream the latest file
# cat "$1" |                        # assumes saved input
#   sed 's/T..:..:..//' |         # hide times
#   cut -d ',' -f 3,1 |           # keep only day and bus ID
#   sort -u |                     # removing duplicate day-buses
#   cut -d ',' -f 2 |             # keep only bus ID
#   sort |                        # preparing for uniq
#   uniq -c |                     # count unique dates
#   sort -k 1 -n |                   # sort in reverse numerical order
#   awk "{print \$2,\$1}"     # print first date, then count

# diff out{1,}

# Using GNU parallel:

INPUT="$1"

process_chunk() {
  local chunk="$1"
  sed 's/T..:..:..//' "$chunk" |
  cut -d ',' -f 3,1                      
}
export -f process_chunk

lines=$(wc -l < "$1")
nproc=$(nproc)
chunk_size=$((lines / nproc))


split -l "$chunk_size" "$INPUT" chunk_
ls chunk_* | parallel -j "$(nproc)" process_chunk > combined.tmp

# Combine and process the results sequentially
cat combined.tmp |
  sort -u |                       # global deduplication
  cut -d ',' -f 2 |              
  sort |                         
  uniq -c |                      
  sort -k 1 -n |                 
  awk '{print $2,$1}'           

rm chunk_*
rm combined.tmp