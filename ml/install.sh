#!/bin/bash

sudo apt-get update

sudo apt-get install -y --no-install-recommends \
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
    joblib==1.4.2 \
    numpy==1.26.4 \
    scikit-learn==1.5.0 \
    scipy==1.13.1 \
    threadpoolctl==3.5.0 \
    imbalanced-learn==0.13.0