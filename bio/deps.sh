#!/bin/bash
# install dependencies
required_version="1.7"

# Check if Samtools is already installed and matches the required version
if command -v samtools &>/dev/null; then
    installed_version=$(samtools --version | head -n 1 | awk '{print $2}')
    if [[ "$installed_version" == "$required_version" ]]; then
        echo "Samtools version $required_version is already installed."
    else
        echo "A different version of Samtools is installed: $installed_version."
        echo "Proceeding to install the required version: $required_version."
    fi
else
    echo "Samtools is not installed. Proceeding with the installation."
    # Update and install prerequisites
    echo "Installing prerequisites..."
    sudo apt update
    sudo apt install -y build-essential libncurses5-dev libncursesw5-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev wget zlib1g-dev

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

    # Verify the installation
    echo "Verifying the installation..."
    installed_version=$(samtools --version | head -n 1 | awk '{print $2}')
    if [[ "$installed_version" == "$required_version" ]]; then
        echo "Samtools version $required_version has been successfully installed."
    else
        echo "Failed to install the correct version of Samtools."
        exit 1
    fi
fi

