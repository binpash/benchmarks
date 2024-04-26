#!/bin/bash

# exit when any command fails
#set -e

IN=$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/
OUT=$PASH_TOP/evaluation/benchmarks/dependency_untangling/output/
IN_NAME=$PASH_TOP/evaluation/benchmarks/dependency_untangling/input/100G.txt

setup_dataset() {
  # download the packages for the package building
  if [ ! -f ${IN}/packages ]; then
      cd $IN
      wget http://pac-n4.csail.mit.edu:81/pash_data/packages
      if [ "$1" = "--small" ]; then
          head -n 20 packages > p
          mv p  packages
      fi
      echo "Package datset downloaded"
  fi
}

source_var() {
  export IN=
}
