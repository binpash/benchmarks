#!/bin/bash

for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done

cd "$(realpath "$(dirname "$0")")" || exit 1
rm -rf ./outputs

if [ "$force" = true ]; then
    rm -rf ./inputs
fi
