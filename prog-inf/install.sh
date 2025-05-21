#!/bin/bash
set -x

TOP=$(git rev-parse --show-toplevel)
URL="https://atlas.cs.brown.edu/data"
installdir="$TOP/prog-inf/inputs"

mkdir -p "$installdir"
cd "$installdir" || exit 1
# Install mir-sa
if [ ! -d mir-sa ]; then
  wget "$URL/prog-inf/mir-sa.tar.gz" -O mir-sa.tar.gz
  tar xf mir-sa.tar.gz
  rm mir-sa.tar.gz
fi

pkgs="node"
for pkg in $pkgs; do
  if ! dpkg -s "$pkg" > /dev/null 2>&1 ; then
    sudo apt-get install -y "$pkg"
  fi
done
