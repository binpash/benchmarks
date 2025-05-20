#!/bin/bash
set -x

TOP=$(git rev-parse --show-toplevel)
URL="https://atlas.cs.brown.edu/data"
installdir="$TOP/prog-inf/inputs"

# Install mir-sa
mkdir -p "$installdir"
cd "$installdir" || exit 1
wget "$URL/prog-inf/mir-sa.zip"
unzip mir-sa.zip
rm mir-sa.zip
