#!/bin/bash

sudo apt update 

sudo apt install -y \
    sudo \
    curl \
    wget \
    unzip \
    python3-pip \
    vim \
    git \
    python3 \
    python3-pip \
    python3-venv

python3 -m venv venv
source venv/bin/activate

pip install --break-system-packages --upgrade pip

pip install --break-system-packages \
    numpy \
    matplotlib