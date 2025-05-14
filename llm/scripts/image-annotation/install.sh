#!/bin/bash

sudo apt-get update

sudo apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    coreutils findutils sed unzip curl imagemagick

pip install --break-system-packages --upgrade pip
pip install --break-system-packages llm
llm install llm-ollama
curl -fsSL https://ollama.com/install.sh | sh
nohup ollama serve > ollama_serve.log 2>&1 &
SERVER_PID=$!
echo "Ollama server started with PID $SERVER_PID."
sleep 5
ollama pull gemma3
