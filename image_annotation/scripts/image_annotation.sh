#!/bin/bash
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/image_annotation"
inputs_dir="$eval_dir/inputs"
outputs_dir="$eval_dir/outputs"

IN=${1:-"$inputs_dir"}
OUT="${2:-"$outputs_dir"}"
mkdir -p "$OUT"

while read -r img; do
    title=$(llm -m gemma3 \
        "Your only output should be a **single** small title for this image:" \
        -a "$img" -o seed 0 < /dev/null)
        
    filename=$(echo "$title" | tr '[:upper:]' '[:lower:]' | sed 's/ /_/g; s/[^a-z0-9_-]//g')
    
    cp "$img" "$OUT/${filename}.jpg"
done < <(find "$IN" -type f -iname "*.jpg" | head -n 10)