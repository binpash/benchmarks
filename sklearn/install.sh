#!/bin/bash

sudo apt-get update

TOP=$(git rev-parse --show-toplevel)
cd "$TOP"/sklearn || exit 1

pip install -r requirements.txt --break-system-packages
