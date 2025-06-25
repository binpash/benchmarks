#!/bin/bash

sudo apt-get update 

pkgs="curl wget unzip coreutils gzip gawk sed findutils git python3 python3-pip python3-venv"

for pkg in $pkgs; do
    if ! dpkg -l | grep -q "$pkg"; then
        sudo apt-get install -y --no-install-recommends "$pkg"
    fi
done

pip install --break-system-packages --upgrade pip

pip install --break-system-packages \
    numpy \
    matplotlib