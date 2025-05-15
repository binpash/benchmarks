#!/bin/bash

input_json="$1"
mrt_file="$2"
output_directory="$3"

zannotate="$(command -v zannotate)"
jq="$(command -v jq)"

annotated="$output_directory/annotated.json"
as_popularity="$output_directory/as_popularity.csv"
file1="$output_directory/ips.txt"
file2="$output_directory/asns.txt"

cat "$input_json" | "$zannotate" -routing -routing-mrt-file="$mrt_file" -input-file-type=json > "$annotated"

cat "$annotated" | "$jq" ".ip" | tr -d '"' > "$file1"
cat "$annotated" | "$jq" -c ".zannotate.routing.asn" > "$file2"

pr -mts, "$file1" "$file2" | awk -F',' '{ a[$2]++; } END { for (n in a) print n "," a[n] }' | sort -k2 -n -t',' -r > "$as_popularity"