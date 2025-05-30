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
