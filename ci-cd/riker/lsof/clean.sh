#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
input_dir="${REPO_TOP}/ci-cd/inputs/scripts/lsof"

rm -rf "$input_dir/dev"
