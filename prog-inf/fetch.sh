#!/bin/bash

TOP=$(git rev-parse --show-toplevel)
URL="https://atlas.cs.brown.edu/data"

size=full
for arg in "$@"; do
    case "$arg" in
    --small) size=full ;; # small uses a subset of full inputs
    --min) size=min ;;
    esac
done

input_dir="${TOP}/prog-inf/inputs"
mkdir -p "$input_dir"
cd "$input_dir" || exit 1

# Fetch the package index
curl -s "${URL}/prog-inf/index.txt" -o index.txt

if [ "$size" = "min" ]; then
  n_packages=2
elif [ "$size" = "small" ]; then
  n_packages=100
else
  n_packages=$(wc -l < index.txt)
fi

# Fetch the packages
mkdir -p "$input_dir"/node

while read -r package; do
  url="$URL/prog-inf/node/$package"
  wget "$url" -O "$input_dir"/node/"$package"
done < index.txt
