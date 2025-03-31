#!/bin/bash

cd "$(dirname "$(realpath "$0")")"

IN="inputs"
IN_NAME="input.txt"

for arg in "$@"; do
    case "$arg" in
        --small) IN_NAME="input_small.txt" ;;
        --min)   IN_NAME="input_min.txt" ;;
    esac
done

if [[ "${1:-}" == "-c" ]]; then
    rm -f ./*.bam ./*.sam
    rm -rf ../output
    exit 0
fi

mkdir -p "$IN" outputs

while IFS= read -r s_line; do
    sample=$(echo "$s_line" | cut -d ' ' -f 2)
    if [[ ! -f "$IN/$sample.bam" ]]; then
        pop=$(echo "$s_line" | cut -d ' ' -f 1)
        link=$(echo "$s_line" | cut -d ' ' -f 3)
        wget -O "$IN/$sample.bam" "$link"
    fi
done < "$IN_NAME"