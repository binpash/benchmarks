#! /bin/sh
OUTPUT=${OUTPUT:-.}
touch "$OUTPUT"/giant
hs_base=$(git rev-parse --show-toplevel)

SCRIPTS="${hs_base}/report/benchmarks/micro/scripts"
"$hs_base/deps/try/try" -y python3 "$SCRIPTS"/giant_file.py "$OUTPUT"/giant 100
