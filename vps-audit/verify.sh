#!/bin/bash

# Exit immediately if a command exits with a non-zero status
# set -e

if [[ "$1" == "--generate" ]]; then
    python3 verify.py --generate
else
    python3 verify.py
fi

echo "vps-audit $?"