#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
IN="$REPO_TOP/aurpkg/input"

mkdir -p "${IN}/deps"

pkgs='ffmpeg unrtf imagemagick libarchive-tools libncurses5-dev libncursesw5-dev zstd liblzma-dev libbz2-dev zip unzip nodejs tcpdump makedeb'

if ! command -v gpg >/dev/null; then
    sudo apt-get install -y gpg
fi

# Add the makedeb repository
wget -qO - 'https://proget.makedeb.org/debian-feeds/makedeb.pub' | gpg --dearmor | sudo tee /usr/share/keyrings/makedeb-archive-keyring.gpg > /dev/null
echo 'deb [signed-by=/usr/share/keyrings/makedeb-archive-keyring.gpg arch=all] https://proget.makedeb.org/ makedeb main' | sudo tee /etc/apt/sources.list.d/makedeb.list > /dev/null
sudo apt update

for pkg in $pkgs; do
    if ! dpkg -s "$pkg" >/dev/null 2>&1; then
        echo "Installing package: $pkg"
        sudo apt-get install -y "$pkg"
    fi
done