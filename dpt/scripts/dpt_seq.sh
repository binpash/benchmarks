#!/bin/bash
# tag: segment-classify-index sequential
# inputs: $1=absolute source image dir, $2=output file path

IMG_DIR="$1"
OUT="$2"

# Limit thread count for determinism
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

pure_func() {
    local img="$1"
    local tmp_boxes
    tmp_boxes=$(mktemp)

    cat "$img" | python3 scripts/sam_segment.py > "$tmp_boxes"
    python3 scripts/classify.py "$img" < "$tmp_boxes" | while read -r line; do
        echo "\"$(basename "$img")\" $line"
    done

    rm "$tmp_boxes"
}
export -f pure_func

rm -f "$OUT"
for img in "$IMG_DIR"/*.png; do
    pure_func "$img" >> "$OUT"
done