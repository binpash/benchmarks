#!/bin/bash

export SUITE_DIR=$(realpath $(dirname "$0"))
cd $SUITE_DIR

# Compression methods to benchmark
COMPRESSIONS=(
  "--gzip"
  "--bzip2"
  "--xz"
  "--lzo"
  "--lz4"
  "--zstd"
)

MAKSELF="./scripts/makeself.sh"

echo "Executing Makeself benchmarks $(date)"

mkdir -p outputs
for compression in "${COMPRESSIONS[@]}"; do
  comp_name=$(echo "$compression" | sed 's/--//')
  input_dir="./inputs/${comp_name}"
  archive_file="./outputs/${comp_name}.run"

  echo "Creating archive with $compression"

  $MAKSELF $compression --current "$input_dir" "$archive_file" \
    "Makeself Archive with $comp_name" ./startup.sh

  if [[ $? -eq 0 ]]; then
    echo "Successfully created $archive_file"
  else
    echo "Failed to create $archive_file"
  fi
done
