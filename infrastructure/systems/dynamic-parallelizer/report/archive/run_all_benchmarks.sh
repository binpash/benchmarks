#!/bin/bash

export PATH=$PATH:$HOME/.local/bin
export PASH_SPEC_TOP=${PASH_SPEC_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

echo "Running all benchmarks"

echo "Running unix 50"
python3 $PASH_SPEC_TOP/report/run_benchmarks.py --subset unix_50 --no-plots
echo "Running transit analytics"
python3 $PASH_SPEC_TOP/report/run_benchmarks.py --subset bus-analytics --no-plots
echo "Running dgsh"
python3 $PASH_SPEC_TOP/report/run_benchmarks.py --subset dgsh --no-plots
echo "Running max temp"
python3 $PASH_SPEC_TOP/report/run_benchmarks.py --subset max_temp --no-plots
