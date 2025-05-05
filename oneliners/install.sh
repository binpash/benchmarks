#!/bin/bash

pkgs="wget bsdmainutils file dos2unix"

sudo apt update

for pkg in $pkgs; do
    if ! dpkg -s "$pkg" >/dev/null 2>&1; then
        echo "Installing package: $pkg"
        sudo apt-get install -y "$pkg"
    fi
done
