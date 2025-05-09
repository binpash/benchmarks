#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# Install Kaggle API client
echo "Installing Kaggle API client..."
pip install kaggle

# Create a .kaggle directory in the home directory
echo "Creating .kaggle directory..."
mkdir -p ~/.kaggle

# Create kaggle.json file with your API credentials
echo "Creating kaggle.json file..."

###  Replace with your own API token  ###
cat > ~/.kaggle/kaggle.json <<EOF
{
  "username": "XXXXXX",
  "key": "XXXXXX"
}
EOF

# Set correct permissions for kaggle.json
chmod 600 ~/.kaggle/kaggle.json

echo "Setup complete!"
