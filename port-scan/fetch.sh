#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
input_dir="${REPO_TOP}/port-scan/inputs"
eval_dir="${REPO_TOP}/port-scan"

size="full"
for arg in "$@"; do
    case "$arg" in
    --small) size="small" ;;
    --min) size="min" ;;
    esac
done

if [ "$size" = "small" ]; then
    N=400000
elif [ "$size" = "min" ]; then
    N=4000
else
    N=40000000
fi
input_dir="$input_dir/$size"
mkdir -p "$input_dir"

if [ -f "$input_dir/port-scan.log" ]; then
    echo "port-scan.log already exists, skipping generation."
else 
    python3 $eval_dir/scripts/generate_inputs.py "$N" > "$input_dir/port-scan.log"
fi
cd "$eval_dir/inputs" || exit 1
if [ -f "routeviews.mrt" ]; then
    echo "routeviews.mrt already exists, skipping download."
    exit 0
fi
wget https://archive.routeviews.org/route-views.linx/bgpdata/2019.10/RIBS/rib.20191001.0000.bz2
bunzip2 rib.20191001.0000.bz2
rm rib.20191001.0000.bz2
mv "${REPO_TOP}/port-scan/inputs/rib.20191001.0000" "${REPO_TOP}/port-scan/inputs/routeviews.mrt"
