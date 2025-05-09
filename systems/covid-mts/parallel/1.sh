#!/bin/bash
# Vehicles on the road per day

# Using GNU parallel:

INPUT="$1"


process_chunk() {
  sed 's/T..:..:..//'|
  cut -d ',' -f 1,3
}
export -f process_chunk

tmp_dir=$(mktemp -d)
trap "rm -rf $tmp_dir" EXIT  

cat "$INPUT" | parallel --pipe --block "$chunk_size" -j "$nproc" process_chunk > "$tmp_dir/combined.tmp"

sort -u "$tmp_dir/combined.tmp" |
  cut -d ',' -f 1 |               
  sort |                          
  uniq -c |                       
  awk '{print $2,$1}'             
