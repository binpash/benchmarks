#!/bin/bash

sudo apt update

sudo apt install -y --no-install-recommends \
    python3.10 \
    python3.10-venv \
    python3.10-distutils \
    python3-pip \
    wget \
    unzip \
    git \
    libgl1 \
    libglib2.0-0 \
    libjpeg-dev \
    zstd \
    ffmpeg \
    imagemagick

python3.10 -m venv .venv
source .venv/bin/activate

pip install --upgrade pip

pip install \
    numpy \
    torch \
    torchvision \
    Pillow \
    segment-anything \
    tensorflow
