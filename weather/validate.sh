#!/bin/bash

set -e

TOP=$(git rev-parse --show-toplevel)

eval_dir="${TOP}/weather"

size="full"

generate=false
for arg in "$@"; do
    case "$arg" in
    --generate) generate=true ;;
    --small) size="small" ;;
    --min) size="min" ;;
    esac
done

statistics_dir="${eval_dir}/outputs/statistics.$size"
correct_dir="${eval_dir}/correct-results/statistics.$size"

if $generate; then
    mkdir -p "$correct_dir"
    cp -r "$statistics_dir"/* "$correct_dir"
    exit 0
fi

diff -q "$statistics_dir/average.txt" "$correct_dir/average.txt"
echo average.$size $?

diff -q "$statistics_dir/min.txt" "$correct_dir/min.txt"
echo min.$size $?

diff -q "$statistics_dir/max.txt" "$correct_dir/max.txt"
echo max.$size $?

hash_dir="$eval_dir/hashes/$size"
hash_file="$hash_dir/tuft-weather.hash"
plot_root="$eval_dir/outputs/$size/plots"

mkdir -p "$hash_dir"

if [ "$generate" = true ]; then
    find "$plot_root" -type f -name '*.png' ! -path '*/tmp/*' -print0 |
        sort -z | tr '\0' '\n' > "$hash_file"
    exit 0
fi

all_exist=true
while IFS= read -r filepath; do
    if [ ! -f "$filepath" ]; then
        echo "Missing: $filepath"
        all_exist=false
    fi
done < "$hash_file"

if [ "$all_exist" = true ]; then
    echo "tuft-weather 0"
else
    echo "tuft-weather 1"
fi
