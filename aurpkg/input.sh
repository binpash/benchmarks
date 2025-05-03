#!/bin/bash

TOP=$(git rev-parse --show-toplevel)
IN="$TOP/aurpkg/inputs"
URL='https://atlas.cs.brown.edu/data'

cd "$TOP" || exit 1
mkdir -p "${IN}"

if [ ! -f "${IN}/packages" ]; then
  wget "$URL"/packages --no-check-certificate -O "${IN}/packages"
fi

for arg in "$@"; do
  if [ "$arg" = "--small" ]; then
    head -n 10 "${IN}/packages" > "${IN}/packages_small"
    break
  elif [ "$arg" = "--min" ]; then
    head -n 1 "${IN}/packages" > "${IN}/packages_min"
    break
  fi
done
