#!/bin/bash
# Hours each vehicle is on the road

# Using GNU parallel:

INPUT="$1"
chunk_size=$(( $(wc -c < "$INPUT") / $(nproc) ))
nproc=$(nproc)
export nproc
export chunk_size
process_chunk() {
  sed 's/T\(..\):..:../,\1/' |  
  cut -d ',' -f 1,2,4                  
}
export -f process_chunk

tmp_dir=$(mktemp -d)
trap "rm -rf $tmp_dir" EXIT  

cat "$INPUT" |
  parallel --pipe --block "$chunk_size" -j "$nproc" process_chunk > "$tmp_dir/combined.tmp"

sort -u "$tmp_dir/combined.tmp" |
  cut -d ',' -f 3 |               
  sort |                         
  uniq -c |                      
  sort -k 1 -n |                 
  awk '{print $2,$1}'             
