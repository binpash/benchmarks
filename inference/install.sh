#!/bin/bash

TOP=$(git rev-parse --show-toplevel)
eval_dir="${TOP}/inference"

sudo apt-get update

sudo apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    libgl1 \
    libglib2.0-0 \
    libjpeg-dev \
    zstd \
    ffmpeg \
    coreutils findutils wget sed unzip curl jq coreutils findutils sed unzip curl imagemagick


pip install --break-system-packages --upgrade pip
pip install --break-system-packages llm
pip install --break-system-packages llm-interpolate
pip install --break-system-packages llm-clap
pip install --break-system-packages llm-ollama

pip install --break-system-packages numpy \
    torch \
    torchvision \
    Pillow \
    segment-anything \
    tensorflow \
    opencv-python

# check if ollama is installed
if ! command -v ollama &> /dev/null
then
    echo "Ollama could not be found, installing..."
    curl -fsSL https://ollama.com/install.sh | sh
else
    echo "Ollama is already installed."
fi
ollama serve > ollama_serve.log 2>&1 &
sleep 5
ollama pull gemma3

ollama_pid=$(pgrep ollama)
kill $ollama_pid