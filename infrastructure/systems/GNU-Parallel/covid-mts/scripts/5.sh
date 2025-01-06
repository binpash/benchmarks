#!/bin/bash
# Hours each bus is active each day

# Records are day, hour, line, bus
# <in.csv sed 's/T\(..\):..:../,\1/' | awk -F, '
# !seen[$1 $2 $4] { seen[$1 $2 $4] = 1; hours[$1 $4]++; bus[$4] = 1; day[$1] = 1; }
# END {
#    PROCINFO["sorted_in"] = "@ind_str_asc"
#    for (d in day)
#      printf("\t%s", d);
#    printf("\n");
#    for (b in bus) {
#      printf("%s", b);
#      for (d in day)
#        printf("\t%s", hours[d b]);
#      printf("\n");
#    }
# }' > out

# Using GNU parallel:

INPUT="$1"

# Function to process chunks
process_chunk() {
  local chunk="$1"
  sed 's/T\(..\):..:../,\1/' "$chunk" | \
  awk -F, '
  !seen[$1 $2 $4] { seen[$1 $2 $4] = 1; print $1, $2, $4; }'
}

export -f process_chunk

# Split the input file into chunks
split -l 10000 "$INPUT" chunk_

# Process each chunk in parallel and combine results
ls chunk_* | parallel -j "$(nproc)" process_chunk > combined.tmp

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
}' combined.tmp > out

# Clean up temporary files
rm chunk_*
rm combined.tmp