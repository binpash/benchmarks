#!/bin/bash
cd "$(realpath "$(dirname "$0")")" || exit 1

REPO_TOP=$(git rev-parse --show-toplevel)
sudo apt-get update -y
sudo apt-get install -y git autoconf automake libtool build-essential cloc time gawk strace lsof python3 python3-pip python3-venv
pip install --break-system-packages -r "$REPO_TOP/infrastructure/requirements.txt"
