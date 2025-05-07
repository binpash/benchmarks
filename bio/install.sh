#!/bin/bash

required_version="1.7"

# Check if Samtools is already installed and matches the required version
if command -v samtools &>/dev/null; then
    installed_version=$(samtools --version | head -n 1 | awk '{print $2}')
    if [[ "$installed_version" == "$required_version" ]]; then
        echo "Samtools version $required_version is already installed."
        exit 0
    fi
fi

sudo apt-get update
pkgs="build-essential libncurses5-dev libncursesw5-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev wget zlib1g-dev"

for pkg in $pkgs; do
    if ! dpkg -s $pkg; then
        sudo apt-get install -y $pkg
    fi
done

# Download Samtools version 1.7
echo "Downloading Samtools version $required_version..."
wget https://github.com/samtools/samtools/releases/download/$required_version/samtools-$required_version.tar.bz2

# Extract the downloaded file
echo "Extracting Samtools..."
tar -xvjf samtools-$required_version.tar.bz2
cd samtools-$required_version

# Compile and install
echo "Compiling and installing Samtools..."
./configure
make
sudo make install

sudo ln -s /usr/local/bin/samtools /usr/bin/samtools

cd ..
rm -rf samtools-$required_version
rm -rf samtools-$required_version.tar.bz2

# Verify the installation
echo "Verifying the installation..."
installed_version=$(samtools --version | head -n 1 | awk '{print $2}')
if [[ "$installed_version" == "$required_version" ]]; then
    echo "Samtools installed version ($installed_version) is $required_version!"
else
    echo "Samtools installed version ($installed_version) is not $required_version!"
    exit 1
fi
