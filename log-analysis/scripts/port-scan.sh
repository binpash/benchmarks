#!/bin/bash

filename="$1"
mrt_file="$2"
annotated="$3"
file1="$4"
file2="$5"
as_popularity="$6"

cat "$filename" |  zannotate -routing -routing-mrt-file=$mrt_file -input-file-type=json > "$annotated"
cat "$annotated" | jq ".ip" | tr -d '"' > "$file1"
cat "$annotated" | jq -c ".zannotate.routing.asn" > "$file2"
pr -mts, $file1 $file2 | awk -F',' "{ a[\$2]++; } END { for (n in a) print n \",\" a[n] } " | sort -k2 -n -t',' -r > "$as_popularity"