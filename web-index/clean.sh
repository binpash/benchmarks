#!/bin/bash

while getopts "f" opt; do
  case $opt in
    f) force=true ;;
  esac
done

rm -r outputs

if [ "$force" = true ]; then
    rm -r inputs
fi
