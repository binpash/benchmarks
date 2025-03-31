#!/bin/bash

cd "$(realpath "$(dirname "$0")")" || exit 1
mkdir -p inputs
cd inputs || exit 1

git clone https://github.com/rozim/ChessData.git

# TODO: Add small inputs and link for download