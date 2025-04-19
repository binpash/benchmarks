#!/bin/bash
# generate-inputs.sh
mkdir -p inputs
INPUT_DIR="${REPO_TOP}/ray-tracing/inputs"
OUT="${INPUT_DIR}/port-scan.log"
rm -f "$OUT"
for i in {1..10000}; do
    printf '{"ip": "192.0.2.%d"}\n' $((RANDOM % 255)) >> "$OUT"
done

cd inputs
wget https://archive.routeviews.org/route-views.linx/bgpdata/2019.10/RIBS/rib.20191001.0000.bz2
bunzip2 rib.20191001.0000.bz2
mv rib.20191001.0000 routeviews.mrt