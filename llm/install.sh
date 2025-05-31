#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/llm"

sudo apt-get update

sudo apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    coreutils findutils sed unzip curl imagemagick jq coreutils findutils sed unzip curl imagemagick


pip install --break-system-packages --upgrade pip
pip install --break-system-packages llm
pip install --break-system-packages llm-interpolate
pip install --break-system-packages llm-clap
pip install --break-system-packages llm-ollama

# check if ollama is installed
if ! command -v ollama &> /dev/null
then
    echo "Ollama could not be found, installing..."
    curl -fsSL https://ollama.com/install.sh | sh
else
    echo "Ollama is already installed."
fi
nohup ollama serve > ollama_serve.log 2>&1 &
SERVER_PID=$!
echo "Ollama server started with PID $SERVER_PID."
sleep 5
ollama pull gemma3
