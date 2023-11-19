#!/bin/bash

echo "Start"
echo " || "
cd inputs
source ./inputs.sh 
echo " || "
source ./run.sh
echo " || "
source ./verify.sh
echo " || "
source ./cleanup.sh
cd ..
echo " || "
echo "We did it!"
