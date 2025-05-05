#!/bin/bash

sudo apt-get update

pkgs='ffmpeg unrtf imagemagick libarchive-tools libncurses5-dev libncursesw5-dev zstd liblzma-dev libbz2-dev zip unzip nodejs tcpdump'

if ! command -v gpg >/dev/null; then
    sudo apt-get install -y gpg
fi

for pkg in $pkgs; do
    if ! dpkg -s "$pkg" >/dev/null 2>&1; then
        sudo apt-get install -y "$pkg"
    fi
done
