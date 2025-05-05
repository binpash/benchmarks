#!/bin/bash

sudo apt-get update 

pkgs="tcpdump curl wget coreutils diffutils gzip bcftools gawk unzip git"

for pkg in $pkgs; do
    if dpkg -s "$pkg" &> /dev/null; then
        sudo apt-get install -y --no-install-recommends "$pkg"
    fi
done
