#!/bin/bash

cd "$(realpath "$(dirname "$0")")" || exit 1

IN="inputs"
IN_NAME="input.txt"

if [[ "${1:-}" == "-c" ]]; then
    rm -f ./*.bam ./*.sam
    rm -rf ../output
    exit 0
fi

for arg in "$@"; do
    case "$arg" in
        --small) IN_NAME="input_small.txt" ;;
        --min)   IN_NAME="input_min.txt" ;;
    esac
done

mkdir -p "${IN}" outputs

if [[ "$IN_NAME" == "input_min.txt" ]]; then
    cp min_inputs/* "$IN/"
    exit 0;
fi

while IFS= read -r s_line; do
    sample=$(echo "$s_line" | cut -d ' ' -f 2)
    if [[ ! -f "$IN/$sample.bam" ]]; then
        pop=$(echo "$s_line" | cut -d ' ' -f 1)
        link=$(echo "$s_line" | cut -d ' ' -f 3)
        wget -O "${IN}/$sample.bam" "$link"
    fi
done < "$IN_NAME"