#!/bin/bash

pkgs="wget bsdmainutils file dos2unix grep findutils mawk"

sudo apt-get update

for pkg in $pkgs; do
    if ! dpkg -s "$pkg" >/dev/null 2>&1; then
        sudo apt-get install -y --no-install-recommends "$pkg"
    fi
done
