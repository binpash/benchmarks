#!/bin/bash

cd "$(realpath "$(dirname "$0")")" || exit 1

while getopts "f" opt; do
  case $opt in
    f) force=true ;;
  esac
done

IN="inputs"
OUT="outputs"

rm -rf "${OUT}"

if [ "$force" = true ]; then
    rm -rf "${IN}/packages"
fi
