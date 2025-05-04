#!/bin/bash
KOALA_SHELL=${KOALA_SHELL:-bash}
base="$(dirname "$0")"
export BENCHMARK_SCRIPT="$base/scripts/main.sh"
export BENCHMARK_CATEGORY="infrastructure/standards/100-files"
$KOALA_SHELL $BENCHMARK_SCRIPT
