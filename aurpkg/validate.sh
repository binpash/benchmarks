#!/bin/bash

# Exit immediately if a command exits with a non-zero status
# set -e

cd "$(realpath "$(dirname "$0")")" || exit 1

[ ! -d "outputs" ] && echo "Directory 'outputs' does not exist" && exit 1

hash_folder="hashes"

size="full"
generate=false
for arg in "$@"; do
    case "$arg" in
        --generate) generate=true ;;
        --small)    size="small"  ;;
        --min)      size="min"    ;;
    esac
done
directory="outputs"

input="inputs/packages"
if [ "$size" = "small" ]; then
    input="inputs/packages_small"
elif [ "$size" = "min" ]; then
    input="inputs/packages_min"
fi

mkdir -p "$hash_folder/$size"

if $generate; then
    while IFS= read -r pkg || [ -n "$pkg" ]; do
        [ -z "$pkg" ] && continue

        file="$directory/$pkg/$pkg.pacscript"

        hash=$(shasum -a 256 "$file" | awk '{print $1}')

        echo "$hash" >"$hash_folder/$size/$pkg.hash"

        echo "$hash_folder/$size/$pkg.hash $hash"
    done <"$input"
    exit 0
fi

while IFS= read -r pkg || [ -n "$pkg" ]; do
    status=0
    [ -z "$pkg" ] && continue

    log="$directory/$pkg/$pkg.pacscript"
    ref="$hash_folder/$size/$pkg.hash"

    if [ ! -f "$log" ] || [ ! -f "$ref" ]; then
        status=1
        continue
    fi

    cur=$(shasum -a 256 "$log" | awk '{print $1}')
    stored=$(cat "$ref")

    [ "$cur" = "$stored" ] || status=1
    echo "$pkg" "$status"
done <"$input"

echo "aurpkg $status"
