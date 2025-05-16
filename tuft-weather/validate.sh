#!/bin/sh

generate=false
size=full

for arg in "$@"; do
    case "$arg" in
        --generate) generate=true ;;
        --small)    size=small    ;;
        --min)      size=min      ;;
    esac
done

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/tuft-weather"

hash_dir="$eval_dir/hashes/$size"
hash_file="$hash_dir/tuft-weather.hash"
plot_root="$eval_dir/outputs/$size/plots"

mkdir -p "$hash_dir"

if [ "$generate" = true ]; then
    find "$plot_root" -type f -name '*.png' ! -path '*/tmp/*' -print0 \
        | sort -z | xargs -0 sha256sum > "$hash_file"
    echo "hash manifest generated"
    exit 0
fi

if ! sha256sum --quiet --check "$hash_file"; then
    echo "tuft-weather 1"
else
    echo "tuft-weather 0"
fi
