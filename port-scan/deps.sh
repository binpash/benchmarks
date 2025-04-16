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
  git

echo "Installing Go $GO_VERSION"
GO_VERSION=1.24.2

# Download and extract to /usr/local
curl -LO https://go.dev/dl/go${GO_VERSION}.linux-amd64.tar.gz
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/port-scan"
go_install_dir="${eval_dir}/go_install"
mkdir -p $go_install_dir
rm -rf /usr/local/go && tar -C $go_install_dir -xzf go${GO_VERSION}.linux-amd64.tar.gz
export PATH=$PATH:/go_install/go/bin
export GOPATH=$HOME/go
export PATH=$PATH:$GOPATH/bin

go install github.com/zmap/zannotate/cmd/zannotate@latest
