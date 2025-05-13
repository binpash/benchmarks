#!/bin/bash

# Exit immediately if a command exits with a non-zero status
# set -e

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/opt-parallel"
outputs_dir="${eval_dir}/outputs"
hash_folder="${eval_dir}/hashes"
[ ! -d "outputs" ] && echo "Directory '$outputs_dir' does not exist" && exit 1
hash_folder="hashes"
generate=false
for arg in "$@"; do
    case "$arg" in
        --small) size=small ;;
        --min)   size=min ;;
        --generate) generate=true ;;
    esac
done

mkdir -p "$hash_folder"

if $generate; then
    for file in $outputs_dir/*.out; do
        filename=$(basename "$file" .out)
        hash=$(shasum -a 256 "$file" | awk '{print $1}')
        echo "$hash" > "$hash_folder/$filename.hash"
        echo "$hash_folder/$filename.hash $hash"
    done
    exit 0
fi

file="$outputs_dir"/opt-parallel_$size.out
filename=$(basename "$file" .out)
hash=$(shasum -a 256 "$file" | awk '{ print $1 }')
echo "$hash" > "outputs/$filename.hash"
diff "$hash_folder/$filename.hash" "outputs/$filename.hash" > /dev/null
match="$?"
echo "$filename $match"
