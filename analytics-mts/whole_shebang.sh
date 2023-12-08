#!/bin/bash

echo "Start"
echo " || "
cd inputs
source ./inputs.sh --full 
echo " || "
source ./run.sh --full
echo " || "
source ./verify.sh
echo " || "
source ./cleanup.sh
cd ..
echo " || "
echo "We did it!"
