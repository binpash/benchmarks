#!/bin/bash

export PASH_SPEC_TOP=${PASH_SPEC_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

benchmark_dir="$PASH_SPEC_TOP/report/benchmarks/sklearn"

cd "$(realpath $(dirname "$0"))"
mkdir -p "$PASH_SPEC_TOP/report/resources/sklearn"
mkdir -p "$PASH_SPEC_TOP/report/output/sklearn"

# Currently just dumped the entire dataset, but ideally we actually download it

pip install -r requirements.txt
