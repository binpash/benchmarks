#!/usr/bin/env bash

set -e
# 7zip
pkgs='p7zip-full curl wget nodejs unzip npm' 
if ! dpkg -s $pkgs >/dev/null 2>&1 ; then
  sudo apt-get install $pkgs -y
  echo 'Packages Installed'
fi

if ! dpkg -s pandoc > /dev/null 2>&1 ; then
  # since pandoc v.2.2.1 does not support arm64, we use v.3.5
  wget https://github.com/jgm/pandoc/releases/download/3.5/pandoc-3.5-1-$(dpkg --print-architecture).deb
  sudo dpkg -i ./pandoc-3.5-1-$(dpkg --print-architecture).deb
  rm ./pandoc-3.5-1-$(dpkg --print-architecture).deb
fi

if ! dpkg -s nodejs > /dev/null 2>&1 ; then
    # node version 18+ does not need external npm
    curl -fsSL https://deb.nodesource.com/setup_18.x | sudo -E bash -
    sudo apt-get install -y nodejs
fi

cd "$(dirname "$0")/scripts" || exit 1
npm install
# Install the npm packages
npm install natural
cd -
