#!/bin/bash

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

pip install --break-system-packages --upgrade pip

pip install --break-system-packages \
    numpy \
    torch \
    torchvision \
    Pillow \
    segment-anything \
    tensorflow \
    opencv-python