#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
IN=$REPO_TOP/aurpkg/input

cd $REPO_TOP || exit 1

mkdir -p ${IN}

# download the packages for the package building
if [ ! -f ${IN}/packages ]; then
  wget https://atlas.cs.brown.edu/data/packages --no-check-certificate -O ${IN}/packages
  echo "Package dataset downloaded"
fi

if [[ "$@" == *"--small"* ]]; then
  head -n 10 ${IN}/packages > ${IN}/packages_small
  mv ${IN}/packages_small ${IN}/packages
fi
