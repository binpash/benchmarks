#!/bin/bash

sudo apt-get update
sudo apt-get install -y --no-install-recommends  gpg \
    wget \
    git \
    unzip \
    zip \
    zstd \
    libncurses5-dev \
    libncursesw5-dev \
    zstd \
    liblzma-dev \
    libbz2-dev \
    zip \
    unzip \
    nodejs \
    libarchive-tools \
    ffmpeg \
    unrtf \
    imagemagick \
    tcpdump \
    cmake \
    build-essential \
    libssl-dev \
    qtcreator qtbase5-dev qt5-qmake gcc libtirpc-dev \
    make \
    libncurses-dev \
    libsm-dev \
    libice-dev \
    libxt-dev \
    libx11-dev \
    libxdmcp-dev \
    libselinux-dev \
    libtool \
    libtool-bin \
    libreadline-dev

wget -qO - 'https://proget.makedeb.org/debian-feeds/makedeb.pub' | gpg --dearmor | sudo tee /usr/share/keyrings/makedeb-archive-keyring.gpg > /dev/null
echo 'deb [signed-by=/usr/share/keyrings/makedeb-archive-keyring.gpg arch=all] https://proget.makedeb.org/ makedeb main' | sudo tee /etc/apt/sources.list.d/makedeb.list > /dev/null

sudo apt-get update
sudo apt-get install -y --no-install-recommends makedeb

TOP=$(git rev-parse --show-toplevel)
URL="https://atlas.cs.brown.edu/data"
installdir="$TOP/pkg/inputs"

mkdir -p "$installdir"
cd "$installdir" || exit 1
# Install mir-sa
if [ ! -d mir-sa ]; then
  wget "$URL/prog-inf/mir-sa.tar.gz" -O mir-sa.tar.gz
  tar xf mir-sa.tar.gz
  rm mir-sa.tar.gz
fi

cd mir-sa/@andromeda/mir-sa || exit 1
if [ ! -d node_modules ]; then
  npm install
fi

# Install Node.js (18.x) and npm via NodeSource
if ! command -v node > /dev/null 2>&1 ; then
  curl -fsSL https://deb.nodesource.com/setup_18.x | sudo -E bash -
  sudo apt-get install -y nodejs
fi

apt install -y default-jdk