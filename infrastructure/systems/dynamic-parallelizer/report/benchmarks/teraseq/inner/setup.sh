#!/bin/sh

export PASH_SPEC_TOP=${PASH_SPEC_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

# Building Docker Container with hs

docker build -t teraseq20-data "$PASH_SPEC_TOP"/report/benchmarks/teraseq/inner
docker build -t hs/teraseq "$PASH_SPEC_TOP" -f "$PASH_SPEC_TOP"/report/benchmarks/teraseq/inner/Dockerfile.hs
