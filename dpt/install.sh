#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/dpt"
rm -rf "$eval_dir/venv" || true
sudo apt update

sudo apt install -y --no-install-recommends \
    wget \
    unzip \
    git \
    libgl1 \
    libglib2.0-0 \
    libjpeg-dev \
    zstd \
    ffmpeg \
    imagemagick \
    parallel \
    python3 \
    python3-pip \
    python3-venv

python3 -m venv "$eval_dir/venv"

pip install --break-system-packages --upgrade pip

pip install --break-system-packages \
    numpy \
    torch \
    torchvision \
    Pillow \
    segment-anything \
    tensorflow \
    opencv-python