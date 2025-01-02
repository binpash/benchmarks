#!/bin/bash

# Exit immediately if a command exits with a non-zero status
# set -e

# run the Python script
python3 verify.py

# check if the script ran successfully
if [ $? -eq 0 ]; then
  echo "verify.py ran successfully."
  exit 0
else
  echo "verify.py encountered an error."
  exit 1
fi