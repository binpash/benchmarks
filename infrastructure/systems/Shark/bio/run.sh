# create bam files with regions
################### 1KG SAMPLES
IN=inputs
IN_NAME=input.txt
OUT=outputs

if [[ "$@" == *"--small"* ]]; then
    IN_NAME=input_small.txt
fi

BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}

"$BENCHMARK_SHELL" ./scripts/bio.sh "$IN" "$IN_NAME" "$OUT"
