#!/bin/bash

IN=${1:-"../inputs"}
OUT="${2:-"../outputs"}"
mkdir -p "$OUT"

while read -r img; do
    title=$(llm -m gemma3 \
        "Your only output should be a **single** small title for this image:" \
        -a "$img" -o seed 0 < /dev/null)
        
    filename=$(echo "$title" | tr '[:upper:]' '[:lower:]' | sed 's/ /_/g; s/[^a-z0-9_-]//g')
    
    cp "$img" "$OUT/${filename}.jpg"
done < <(find "$IN" -type f -iname "*.jpg")