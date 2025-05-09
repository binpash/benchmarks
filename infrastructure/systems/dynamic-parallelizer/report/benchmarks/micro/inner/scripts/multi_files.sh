#! /bin/sh
OUTPUT=${OUTPUT:-.}
SCRIPTS=${SCRIPTS:-./scripts}
python3 "$SCRIPTS"/multi_files.py "$OUTPUT"/foo
