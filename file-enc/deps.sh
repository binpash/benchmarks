#!/bin/bash

sudo apt-get update

required_pkgs=(
  ffmpeg
  unrtf
  imagemagick
  libarchive-tools
  libncurses5-dev
  libncursesw5-dev
  zstd
  liblzma-dev
  libbz2-dev
  zip
  unzip
  nodejs
  tcpdump
)

missing_pkgs=()

for pkg in "${required_pkgs[@]}"; do
  if ! dpkg -s "$pkg" >/dev/null 2>&1; then
    missing_pkgs+=("$pkg")
  fi
done

if [ "${#missing_pkgs[@]}" -gt 0 ]; then
  echo "Installing missing packages: ${missing_pkgs[*]}"
  sudo apt-get install -y "${missing_pkgs[@]}"
fi