#!/bin/bash

# Directory containing the web-index structure
BENCH_TOP=${BENCH_TOP:-$(git rev-parse --show-toplevel)}
RESOURCE_DIR=${RESOURCES_DIR:-$BENCH_TOP/report/resources/web-index/}

# Inflation command
INFLATE_SCRIPT="python3 $BENCH_TOP/report/util/inflate.py"

rm -rf $RESOURCE_DIR/articles100m500k
mkdir -p $RESOURCE_DIR/articles100m500k
cp -r "$RESOURCE_DIR/articles100m" "$RESOURCE_DIR/articles100m500k"
# Traverse all HTML files and inflate them to 500KB
find "$RESOURCE_DIR/articles100m500k" -type f -name "*.html" | while read -r html_file; do
    cd "$(dirname $html_file)"
    # Inflate the HTML file to 500KB
    $INFLATE_SCRIPT "$html_file" 500K
    filename=$(basename -- "$html_file")
    mv "500K-$filename" "$html_file"
done

rm -rf $RESOURCE_DIR/articles100m1m
mkdir -p $RESOURCE_DIR/articles1g1m
cp -r "$RESOURCE_DIR/articles100m" "$RESOURCE_DIR/articles100m1m"
# Traverse all HTML files and inflate them to 1MB
find "$RESOURCE_DIR/articles100m1m" -type f -name "*.html" | while read -r html_file; do
    cd "$(dirname $html_file)"
    # Inflate the HTML file to 1MB
    $INFLATE_SCRIPT "$html_file" 1M
    filename=$(basename -- "$html_file")
    mv "1M-$filename" "$html_file"
done

# mkdir -p $RESOURCE_DIR/articles1g500k
# cp -r "$RESOURCE_DIR/articles100m" "$RESOURCE_DIR/articles1g500k"
# # Traverse all HTML files and inflate them to 500KB
# find "$RESOURCE_DIR/articles1g500k" -type f -name "*.html" | while read -r html_file; do
#     cd "$(dirname $html_file)"
#     # Inflate the HTML file to 500KB
#     $INFLATE_SCRIPT "$html_file" 500K
#     rm $html_file
# done

# mkdir -p $RESOURCE_DIR/articles1g1m
# cp -r "$RESOURCE_DIR/articles100m" "$RESOURCE_DIR/articles1g1m"
# # Traverse all HTML files and inflate them to 1MB
# find "$RESOURCE_DIR/articles1g1m" -type f -name "*.html" | while read -r html_file; do
#     cd "$(dirname $html_file)"
#     # Inflate the HTML file to 1MB
#     $INFLATE_SCRIPT "$html_file" 1M
#     rm $html_file
# done
