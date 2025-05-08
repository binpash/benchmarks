#!/bin/bash

for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done

rm -r outputs

if [ "$force" = true ]; then
    rm -r inputs
fi
