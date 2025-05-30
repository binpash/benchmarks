#!/bin/bash

GENERATE=false

for arg in "$@"; do
    if [[ "$arg" == "--generate" ]]; then
        GENERATE=true
        break
    fi
done

if [[ "$GENERATE" == true ]]; then
    python3 validate.py --generate
    exit 0
else
    python3 validate.py
    echo "vps-audit $?"
fi
