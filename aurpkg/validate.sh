#!/bin/bash

# Exit immediately if a command exits with a non-zero status
# set -e

cd "$(realpath "$(dirname "$0")")" || exit 1

[ ! -d "outputs" ] && echo "Directory 'outputs' does not exist" && exit 1

size="full"
for arg in "$@"; do
    case "$arg" in
        --small)    size="small"  ;;
        --min)      size="min"    ;;
    esac
done

directory="outputs"
input="inputs/packages"
N=150

if [ "$size" = "small" ]; then
    input="inputs/packages_small"
    N=1
elif [ "$size" = "min" ]; then
    input="inputs/packages_min"
    N=10
fi


missing=0

while IFS= read -r pkg || [ -n "$pkg" ]; do
    file="$directory/$pkg.txt"
    if [ ! -f "$file" ]; then
        missing=$((missing + 1))
        continue
    fi
    if ! grep -q "Finished making" "$file"; then
        missing=$((missing + 1))
    fi
done < "$input"

if [ "$missing" -eq 0 ]; then
    echo "aurpkg 0"
else
    echo "aurpkg 1"
fi