#!/bin/bash

for arg in "$@"; do
    case "$arg" in
        "-f") force=true ;;
    esac
done

rm -rf outputs

if [ "$force" = true ]; then
    rm -rf inputs
fi
