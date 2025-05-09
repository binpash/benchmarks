#!/bin/bash
# Hours monitored each day

# Using GNU parallel:

INPUT="$1"
chunk_size=$(( $(wc -c < "$INPUT") / $(nproc) ))
nproc=$(nproc)
export nproc
export chunk_size
process_chunk() {
  local chunk="$1"
  sed 's/T\(..\):..:../,\1/'|
  cut -d ',' -f 1,2                    
}
export -f process_chunk

tmp_dir=$(mktemp -d)
trap "rm -rf $tmp_dir" EXIT  

cat "$INPUT" |
  parallel --pipe --block "$chunk_size" -j "$nproc" process_chunk > "$tmp_dir/combined.tmp"

sort -u "$tmp_dir/combined.tmp" |
  cut -d ',' -f 1 |               
  sort |                         
  uniq -c |                      
  awk '{print $2,$1}'            
