#!/bin/bash
# Hours each bus is active each day

INPUT="$1"
chunk_size=$(( $(wc -c < "$INPUT") / $(nproc) ))
nproc=$(nproc)
export nproc
export chunk_size
# Function to process chunks
process_chunk() {
  local chunk="$1"
  sed 's/T\(..\):..:../,\1/' "$chunk" | \
  awk -F, '
  !seen[$1 $2 $4] { seen[$1 $2 $4] = 1; print $1, $2, $4; }'
}

export -f process_chunk

tmp_dir=$(mktemp -d)
trap "rm -rf $tmp_dir" EXIT  

# Process each chunk in parallel and combine results
cat "$INPUT" | parallel --pipe --block "$chunk_size" -j "$(nproc)" process_chunk > "$tmp_dir/combined.tmp"

# Aggregate results globally
awk '
{
  if (!seen[$1 " " $2 " " $3]++) {
    hours[$1 " " $3]++; 
    bus[$3] = 1; 
    day[$1] = 1;
  }
}
END {
  PROCINFO["sorted_in"] = "@ind_str_asc"
  # Print header row with days
  printf("\t");
  for (d in day)
    printf("%s\t", d);
  printf("\n");

  # Print bus rows with activity counts
  for (b in bus) {
    printf("%s\t", b);
    for (d in day)
      printf("%d\t", hours[d " " b] ? hours[d " " b] : 0);
    printf("\n");
  }
}' "$tmp_dir/combined.tmp" > out
  