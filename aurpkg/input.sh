#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
IN="$REPO_TOP/aurpkg/inputs"

cd "$REPO_TOP" || exit 1
mkdir -p "${IN}"

if [ ! -f "${IN}/packages" ]; then
  wget https://atlas.cs.brown.edu/data/packages --no-check-certificate -O "${IN}/packages"
fi

for arg in "$@"; do
  if [ "$arg" = "--small" ]; then
    head -n 10 "${IN}/packages" > "${IN}/packages_small"
    mv "${IN}/packages_small" "${IN}/packages"
    break
  elif [ "$arg" = "--min" ]; then
    head -n 1 "${IN}/packages" > "${IN}/packages_min"
    mv "${IN}/packages_min" "${IN}/packages"
    break
  fi
done
