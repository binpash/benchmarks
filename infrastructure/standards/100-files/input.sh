#!/bin/bash

inputs_dir="$(dirname "$0")/inputs"

mkdir -p "$inputs_dir"

for i in {0..99}; do
    echo $i > "$inputs_dir/$(printf "%03d" $i).txt"
done
