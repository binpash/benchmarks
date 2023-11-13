#!/bin/bash

echo "Start"
echo " || "
cd inputs
#source ./inputs.sh
echo " || "
source ./run.sh --small
echo " || "
source ./verify.sh --small
echo " || "
#source ./cleanup.sh
cd ..
echo " || "
echo "We did it!"
