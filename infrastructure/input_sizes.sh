#!/bin/bash
set -euo pipefail

REPO_TOP="$(git rev-parse --show-toplevel)"
cd "$REPO_TOP" || exit 1
SIZE_FLAG="--min"
SUFFIX="min"

KEEP=false
for arg in "$@"; do
    case "$arg" in
    -k) KEEP=true ;;
    --min)
        SIZE_FLAG="--min"
        SUFFIX="min"
        ;;
    --small)
        SIZE_FLAG="--small"
        SUFFIX="small"
        ;;
    --full)
        SIZE_FLAG="--full"
        SUFFIX="full"
        ;;
    *)
        echo "Unknown argument: $arg"
        exit 1
        ;;
    esac
done

CSV_OUT="$REPO_TOP/input_sizes.${SUFFIX}.csv"

default_benchmarks=(
    "analytics"
    "bio"    
    "covid"
    "file-mod"
    "inference"
    "ml"
    "nlp"
    "oneliners"
    "pkg"
    "unixfun"
    "weather"
    "web-search"
    "web-search.new"
)

echo "benchmark,size_bytes" >"$CSV_OUT"

if [[ "$KEEP" == true ]]; then
    echo "Keeping existing inputs"
else
    echo "Purging inputs"
    for bench in "${default_benchmarks[@]}"; do
        dir="$REPO_TOP/$bench"
        if [[ -d "$dir/inputs" ]]; then
            echo "Purging inputs in $bench"
            rm -rf "$dir/inputs"
        fi
    done
fi 

for bench in "${default_benchmarks[@]}"; do
    dir="$REPO_TOP/$bench"
    echo ">>> $bench"

    if [[ ! -d "$dir" ]]; then
        echo "Directory $dir does not exist. Skipping."
        continue
    fi

    cd "$dir" || continue

    if [[ ! -x "./fetch.sh" ]]; then
        echo "No fetch.sh in $bench. Skipping."
        continue
    fi
    ./install.sh
    ./fetch.sh "$SIZE_FLAG" || {
        echo "Fetch failed in $bench"
        continue
    }

    if [[ -d inputs ]]; then
        size_bytes=$(du -sb inputs | awk '{print $1}')
        echo "$bench,$size_bytes" >>"$CSV_OUT"
        rm -rf inputs
    else
        echo "No inputs/ dir in $bench"
    fi
done

echo "CSV saved to: $CSV_OUT"
