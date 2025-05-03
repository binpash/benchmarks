#!/bin/bash

GENERATE=false

for arg in "$@"; do
    if [[ "$arg" == "--generate" ]]; then
        GENERATE=true
        break
    fi
done

if [[ "$GENERATE" == true ]]; then
    python3 verify.py --generate
else
    python3 verify.py
fi

echo "vps-audit $?"