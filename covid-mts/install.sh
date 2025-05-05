#!/bin/bash

sudo apt-get update 

pkgs="coreutils curl gzip gawk sed git"

for pkg in $pkgs; do
    if ! dpkg -s "$pkg" &> /dev/null; then
        sudo apt install -y "$pkg"
    fi
done
