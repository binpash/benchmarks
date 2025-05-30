#! /bin/bash

sudo apt-get update 

pkgs="binutils git build-essential coreutils wget unzip make pbzip2 binutils bzip2 zstd gnupg"

for pkg in $pkgs; do
    if ! dpkg -s "$pkg" &> /dev/null; then
        sudo apt-get install -y "$pkg"
    fi
done
