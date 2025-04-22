#!/bin/bash
BENCHMARK_SHELL=${BENCHMARK_SHELL:-bash}
base="$(dirname "$0")"
export BENCHMARK_SCRIPT="$base/scripts/main.sh"
export BENCHMARK_CATEGORY="infrastructure/standards/time-in-shell-subprocess"
$BENCHMARK_SHELL $BENCHMARK_SCRIPT
