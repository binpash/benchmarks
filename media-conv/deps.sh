#!/bin/bash

sudo apt update 

sudo apt install -y --no-install-recommends \
    sudo \
    curl \
    wget \
    unzip \
    python3-pip \
    vim \
    ffmpeg unrtf imagemagick libarchive-tools libncurses5-dev libncursesw5-dev zstd liblzma-dev libbz2-dev zip unzip nodejs tcpdump \
    git

wget https://download.imagemagick.org/archive/releases/ImageMagick-6.9.11-60.tar.xz
tar -xf ImageMagick-6.9.11-60.tar.xz 
cd ImageMagick-6.9.11-60 || exit 1
./configure --prefix=/opt/imagemagick-6.9.11-60
make -j$(nproc)
sudo make install
export PATH=/opt/imagemagick-6.9.11-60/bin:$PATH