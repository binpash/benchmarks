#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/llm/scripts/playlist-creation"
rm -rf "$eval_dir/venv" || true

sudo apt-get update

sudo apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    coreutils findutils sed unzip curl imagemagick jq

if [ ! -d "$eval_dir/venv" ]; then
    python3 -m venv "$eval_dir/venv"
fi

pip install --break-system-packages --upgrade pip
pip install --break-system-packages llm
pip install --break-system-packages llm-interpolate
pip install --break-system-packages llm-clap