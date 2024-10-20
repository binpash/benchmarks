#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/unix50"
input_dir="${eval_dir}/input"

mkdir -p $input_dir/small
mkdir -p $input_dir/full

for input in 1 10 11 12 2 3 4 5 6 7 8 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9
do
    content="$(curl --insecure "https://atlas-group.cs.brown.edu/data/unix50/${input}.txt")"

    small="$input_dir/small/${input}_1M.txt"
    yes "$content" | head -c 1048576 > $small

    filename="$input_dir/full/${input}_3G.txt"
    truncate -s0 $filename
    for i in {0..1000}; do # can change this back to 3000 for 3G
        cat $small >> $filename
    done
done
