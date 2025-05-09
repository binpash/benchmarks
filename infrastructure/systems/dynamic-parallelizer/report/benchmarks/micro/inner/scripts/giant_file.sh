#! /bin/sh
hs_base=$(git rev-parse --show-toplevel)
OUTPUT=${OUTPUT:-.}
SCRIPTS="${hs_base}/report/benchmarks/micro/inner/scripts"
touch "$OUTPUT"/giant
python3 "$SCRIPTS"/giant_file.py "$OUTPUT"/giant 100
