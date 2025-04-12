#!/bin/bash
set -e

REPO_TOP=$(git rev-parse --show-toplevel)
INPUT_DIR="${REPO_TOP}/ray-tracing/inputs"
mkdir -p "$INPUT_DIR"

# TODO add inputs
# Simulate input fetch (replace with actual source)
# curl -o "$INPUT_DIR/1.INFO" <URL>
# curl -o "$INPUT_DIR/2.INFO" <URL>
echo "Expecting logs manually placed in $INPUT_DIR"
python3 get_inputs.py