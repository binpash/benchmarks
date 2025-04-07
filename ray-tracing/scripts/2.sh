#!/bin/bash

in_dir=$1
out_dir=$2

cat "$in_dir"/*.INFO | grep "\[RAY\]" | grep -v pathID | cut -c 7- >> "$out_dir/rays.csv"
