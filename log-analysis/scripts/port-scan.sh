#!/bin/bash

mrt_file="$2"
output_directory="$3"
mkdir -p "$output_directory"

zannotate="$(command -v zannotate)"
jq="$(command -v jq)"

pure_func() {
    log_file="$1"
    base_name=$(basename "$log_file")
    log_output_dir="$output_directory/$base_name"
    mkdir -p "$log_output_dir"

    zannotate="$(command -v zannotate)"
    jq="$(command -v jq)"

    annotated="$log_output_dir/annotated.json"
    as_popularity="$log_output_dir/as_popularity.csv"
    file1="$log_output_dir/ips.txt"
    file2="$log_output_dir/asns.txt"

    # annotate the log
    cat "$log_file" | "$zannotate" -routing -routing-mrt-file="$mrt_file" -input-file-type=json > "$annotated"

    # extract IPs and ASNs
    cat "$annotated" | "$jq" ".ip" | tr -d '"' > "$file1"
    cat "$annotated" | "$jq" -c ".zannotate.routing.asn" > "$file2"

    # compute ASN popularity
    pr -mts, "$file1" "$file2" | awk -F',' '{ a[$2]++; } END { for (n in a) print n "," a[n] }' | sort -k2 -n -t',' -r > "$as_popularity"
}
export -f pure_func

# process each log
for log in "$1"/*; do
    out_dir="$output_directory/$(basename "$log")"
    mkdir -p "$out_dir"
    pure_func "$log"
done
