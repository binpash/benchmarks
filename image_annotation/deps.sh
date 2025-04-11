#!/bin/bash

sudo apt-get update

sudo apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    coreutils findutils sed unzip curl imagemagick

python3 -m venv venv
source venv/bin/activate

pip install --upgrade pip
pip install llm
llm install llm-ollama
curl -fsSL https://ollama.com/install.sh | sh
ollama serve || true &
ollama pull gemma3