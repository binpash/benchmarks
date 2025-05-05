#!/bin/bash

benchmarks="oneliners bio vps-audit file-enc log-analysis web-index media-conv covid-mts sklearn nlp aurpkg unix50 makeself riker max-temp"

echo "[*] Running kick tires script"
echo "[*] Running ./main.sh --min for each benchmark"

# STEP 1: Run ./main.sh --min in each benchmark directory
for bench in $benchmarks; do
  echo "Running benchmark: $bench"
  ./main.sh --min "$bench" || {
    echo "Failed to run --min for $bench"
    exit 1
  }
done
