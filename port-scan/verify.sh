#!/bin/bash

cd "$(realpath "$(dirname "$0")")" || exit 1

outputs_dir="outputs"
hash_folder="hashes"
generate=false

[ ! -d "$outputs_dir" ] && echo "Directory '$outputs_dir' does not exist" && exit 1

for arg in "$@"; do
    if [[ "$arg" == "--generate" ]]; then
        generate=true
        continue
    fi
    case "$arg" in
        --small) hash_folder="hashes/small" ;;
        --min) hash_folder="hashes/min" ;;
    esac
done

mkdir -p "$hash_folder"

as_popularity_file="$outputs_dir/as_popularity.csv"
annotated_file="$outputs_dir/annotated.json"

if $generate; then
    for file in "$as_popularity_file" "$annotated_file"; do
        [ -f "$file" ] || continue
        filename=$(basename "$file")
        hash=$(shasum -a 256 "$file" | awk '{ print $1 }')
        echo "$hash" > "$hash_folder/$filename.hash"
        echo "$hash_folder/$filename.hash $hash"
    done
    exit 0
fi

all_ok=0
for file in "$as_popularity_file" "$annotated_file"; do
    filename=$(basename "$file")
    hash_file="$hash_folder/$filename.hash"
    actual_hash=$(shasum -a 256 "$file" | awk '{ print $1 }')

    if [[ ! -f "$hash_file" ]]; then
        echo "Missing hash file: $hash_file"
        all_ok=1
        continue
    fi

    expected_hash=$(cat "$hash_file")
    if [[ "$actual_hash" != "$expected_hash" ]]; then
        echo "$filename 1"
        all_ok=1
    else
        echo "$filename 0"
    fi
done

exit $all_ok