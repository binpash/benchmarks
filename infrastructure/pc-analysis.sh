#!/bin/bash

TOP=$(git rev-parse --show-toplevel)
cd "$TOP"/infrastructure || exit 1
OUTPUT="$1"

[ -z "$OUTPUT" ] && {
    echo "Usage: $0 <output_file>"
    exit 1
}

python3 ./make_pca_plot.py "$OUTPUT"
