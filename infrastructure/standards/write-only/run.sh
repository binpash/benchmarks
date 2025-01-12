#!/bin/bash
BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
base="$(dirname "$0")"
export BENCHMARK_SCRIPT="$base/scripts/main.sh"
export BENCHMARK_CATEGORY="infrastructure/standards/write-only"
$BENCHMARK_SHELL $BENCHMARK_SCRIPT
