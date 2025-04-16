#!/bin/bash
set -e

REPO_TOP=$(git rev-parse --show-toplevel)
INPUT_DIR="${REPO_TOP}/ray-tracing/inputs"
eval_dir="${REPO_TOP}/port-scan"
mkdir -p "$INPUT_DIR"

# TODO add inputs
# Simulate input fetch (replace with actual source)
# curl -o "$INPUT_DIR/1.INFO" <URL>
# curl -o "$INPUT_DIR/2.INFO" <URL>
# bash "$REPO_TOP/ray-tracing/generate_inputs.sh" "$INPUT_DIR"

cp $eval_dir/temp_inputs/* "$INPUT_DIR"

cd "$INPUT_DIR" || exit 1
wget https://archive.routeviews.org/route-views.linx/bgpdata/2019.10/RIBS/rib.20191001.0000.bz2
bunzip2 rib.20191001.0000.bz2
mv rib.20191001.0000 routeviews.mrt