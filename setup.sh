#!/bin/bash
set -e

cd "$(realpath "$(dirname "$0")")" || exit 1

TOP=$(git rev-parse --show-toplevel)
sudo apt-get update
sudo apt-get install -y  git autoconf automake libtool build-essential cloc time gawk jq strace lsof python3 python3-pip python3-venv
VENV_DIR="$TOP/venv"
if [ ! -d "$VENV_DIR" ]; then
    python3 -m venv "$VENV_DIR"
fi
source "$VENV_DIR/bin/activate"
pip install --upgrade pip
pip install --break-system-packages -r "$TOP/infrastructure/requirements.txt"
