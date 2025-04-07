#!/bin/bash

in_dir=$1
out_dir=$2

cat "$in_dir"/*.INFO | grep "\[RAY\]" | grep -v pathID | cut -c 7- | sed -n '/^590432,/p' >> "$out_dir/rays.csv"
