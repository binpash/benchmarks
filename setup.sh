#!/bin/bash
cd "$(realpath "$(dirname "$0")")" || exit 1

TOP=$(git rev-parse --show-toplevel)
sudo apt-get update
sudo apt-get install -y  git autoconf automake libtool build-essential cloc time gawk jq strace lsof python3 python3-pip python3-venv
pip install --break-system-packages -r "$TOP/infrastructure/requirements.txt"
