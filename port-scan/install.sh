#!/bin/bash

sudo apt-get update 

sudo apt-get install -y --no-install-recommends \
  jq \
  coreutils \
  gawk \
  cmake \
  build-essential \
  libjansson-dev \
  libpcap-dev \
  tar \
  git \
  python3

echo "Installing Go $GO_VERSION"
GO_VERSION=1.24.2

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/port-scan"
go_install_dir="${eval_dir}/go_install"

mkdir -p "$go_install_dir"
curl -LO "https://go.dev/dl/go${GO_VERSION}.linux-amd64.tar.gz"
tar -C "$go_install_dir" -xzf "go${GO_VERSION}.linux-amd64.tar.gz"
rm -f "go${GO_VERSION}.linux-amd64.tar.gz"

export GOROOT="$go_install_dir/go"
export GOPATH="$HOME/go"
export PATH="$GOROOT/bin:$GOPATH/bin:$PATH"

go install github.com/zmap/zannotate/cmd/zannotate@latest
