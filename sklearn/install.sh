#!/bin/bash

sudo apt update

export PASH_SPEC_TOP=${PASH_SPEC_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

cd "$(realpath "$(dirname "$0")")" || exit 1
mkdir -p "$PASH_SPEC_TOP/report/resources/sklearn"
mkdir -p "$PASH_SPEC_TOP/report/output/sklearn"

# Currently just dumped the entire dataset, but ideally we actually download it

pip install -r requirements.txt --break-system-packages
