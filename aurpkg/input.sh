#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
IN=$REPO_TOP/aurpkg/input

# download the packages for the package building
if [ ! -f ${IN}/packages ]; then
  cd $IN
  wget https://atlas.cs.brown.edu/data/packages --no-check-certificate
  if [ "$1" = "--small" ]; then
      head -n 20 packages > p
      mv p  packages
  fi
  echo "Package datset downloaded"
fi
