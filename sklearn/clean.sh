#!/bin/bash

cd "$(realpath "$(dirname "$0")")" || exit 1

for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done

rm -rf ./result

if [ "$force" = true ]; then
    rm -rf ./tmp
    rm -rf ./inputs
fi
