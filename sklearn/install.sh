#!/bin/bash

sudo apt-get update

TOP=$(git rev-parse --show-toplevel)
cd "$TOP"/sklearn || exit 1

# check if recent-enough pip version
if pip install --help | grep -q -- '--break-system-packages'; then
    true
else
    # upgrade pip
    python3 -m pip install --upgrade pip --user >/dev/null 2>&1
fi

pip install -r requirements.txt --break-system-packages
