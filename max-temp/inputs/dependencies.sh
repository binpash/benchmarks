#!/usr/bin/env bash

echo "Clearing depedencies"

pkgs='curl wget' 
if ! dpkg -s $pkgs >/dev/null 2>&1 ; then
  sudo apt-get install $pkgs -y
  echo 'Packages Installed'
fi
