#!/bin/bash

sudo apt-get update

pkgs='ffmpeg unrtf imagemagick libarchive-tools libncurses5-dev libncursesw5-dev zstd liblzma-dev libbz2-dev zip unzip nodejs tcpdump'

if ! dpkg -s "$pkgs" >/dev/null 2>&1 ; then
    sudo apt-get install "$pkgs" -y
fi
