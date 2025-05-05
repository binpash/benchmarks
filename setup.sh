#!/bin/bash
cd "$(realpath "$(dirname "$0")")" || exit 1

sudo apt-get install -y autoconf automake libtool build-essential cloc time gawk python3 python3-pip python3-venv
python3 -m venv venv
source venv/bin/activate
pip install --break-system-packages -r "$REPO_TOP/infrastructure/requirements.txt"
