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


REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/media-conv"
install_dir="${eval_dir}/install"
mkdir -p "$install_dir"
cd "$install_dir" || exit 1
wget https://download.imagemagick.org/archive/releases/ImageMagick-6.9.11-60.tar.xz
tar -xf ImageMagick-6.9.11-60.tar.xz 
cd ImageMagick-6.9.11-60 || exit 1
./configure --prefix="$install_dir/imagemagick-6.9.11-60/install"
make -j$(nproc)
sudo make install
export PATH="$install_dir/ImageMagick-6.9.11-60/install/bin:$PATH"
