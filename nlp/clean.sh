#!/bin/bash

while getopts "f" opt; do
  case $opt in
    f) force=true ;;
  esac
done

cd "$(realpath "$(dirname "$0")")" || exit 1

rm -rf ./outputs

if [ "$force" = true ]; then
    rm -rf ./inputs
fi
