#!/bin/bash

in_dir=$1
out_dir=$2

cat "$in_dir/1.INFO" | grep "\[RAY\]" | head -n1 | cut -c 7- > "$out_dir/rays.csv"
