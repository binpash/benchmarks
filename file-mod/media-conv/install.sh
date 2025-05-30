#!/bin/bash

sudo apt-get update 

pkgs="curl wget unzip gzip coreutils ffmpeg unrtf imagemagick zstd git xz-utils"

for pkg in $pkgs; do
    if ! dpkg -s "$pkg" &> /dev/null; then
        sudo apt-get install -y --no-install-recommends "$pkg"
    fi
done
