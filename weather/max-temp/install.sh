#!/bin/bash

sudo apt-get update 

pkgs="curl wget coreutils gzip gawk sed findutils git"

for pkg in $pkgs; do
    if ! dpkg -l | grep -q "$pkg"; then
        sudo apt-get install -y "$pkg"
    fi
done
