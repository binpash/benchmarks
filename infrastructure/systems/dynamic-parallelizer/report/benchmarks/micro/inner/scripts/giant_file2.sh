#! /bin/sh
OUTPUT=${OUTPUT:-.}
hs_base=$(git rev-parse --show-toplevel)
SCRIPTS="${hs_base}/report/benchmarks/micro/scripts"
touch "$OUTPUT"/giant
python3 "$SCRIPTS"/giant_file.py "$OUTPUT"/giant 10000
