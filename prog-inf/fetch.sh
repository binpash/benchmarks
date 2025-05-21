#!/bin/bash
set -x
set -e

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
curl -s "${URL}/prog-inf/node_modules/index.txt" -o index.txt

if [ "$size" = "min" ]; then
  n_packages=2
elif [ "$size" = "small" ]; then
  n_packages=100
else
  n_packages=$(wc -l < index.txt)
fi
head -n "$n_packages" <index.txt > index.$size.txt

# Fetch the packages
mkdir -p "$input_dir"
if [ ! -d node_modules ]; then
  wget "$URL/prog-inf/node_modules.tar.gz" -O node_modules.tar.gz
  tar xf node_modules.tar.gz
  rm node_modules.tar.gz
fi

