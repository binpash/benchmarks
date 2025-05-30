#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir=${REPO_TOP}/tuft-weather
sudo apt update

sudo apt install -y --no-install-recommends \
    wget \
    unzip \
    gawk \
    python3 \
    python3-pip \
    python3-venv

pip install --break-system-packages --upgrade pip

pip install --break-system-packages \
    numpy \
    matplotlib