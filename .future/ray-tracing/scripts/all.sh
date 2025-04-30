#!/bin/bash
# source: posh benchmark suite

in_dir=$1
out_dir=$2

mkdir -p "$out_dir"

cat "$in_dir/1.INFO" | grep "\[RAY\]" | head -n1 | cut -c 7- > "$out_dir/rays.csv"

cat "$in_dir"/*.INFO | grep "\[RAY\]" | grep -v pathID | cut -c 7- >> "$out_dir/rays.csv"

cat "$in_dir"/*.INFO | grep "\[RAY\]" | grep -v pathID | cut -c 7- | sed -n '/^590432,/p' >> "$out_dir/rays.csv"

cat "$out_dir/rays.csv" | q -H -d, "SELECT * FROM - WHERE pathID = 20613314" > "$out_dir/query_pathid.log"

cat "$out_dir/rays.csv" | q -H -d, "SELECT MAX(timestamp), MAX(hop) FROM - GROUP BY pathID LIMIT 5" > "$out_dir/query_max_hop.log"
