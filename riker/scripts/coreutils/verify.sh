#!/bin/sh
REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input"
COREUTILS_DEV="$input_dir/scripts/coreutils/dev"

COREUTILS_DIR="$COREUTILS_DEV/src"

# List of core utilities to test
TOOLS="echo cat ls head tail"

# Verify tools exist
for tool in $TOOLS; do
  TOOL_PATH="$COREUTILS_DIR/$tool"
  if [ ! -x "$TOOL_PATH" ]; then
    echo riker/coreutils 1; 
    exit 1
  fi
done

# Run functional checks

# echo
OUT=$("$COREUTILS_DIR/echo" "hello")
[ "$OUT" = "hello" ] || { echo riker/coreutils 1; exit 1; }

# cat
echo "test" > tmpfile.txt
OUT=$("$COREUTILS_DIR/cat" tmpfile.txt)
[ "$OUT" = "test" ] || { rm -f tmpfile.txt; echo riker/coreutils 1; exit 1; }

# head
OUT=$("$COREUTILS_DIR/head" -n 1 tmpfile.txt)
[ "$OUT" = "test" ] || { rm -f tmpfile.txt; echo riker/coreutils 1; exit 1; }

# tail
OUT=$("$COREUTILS_DIR/tail" -n 1 tmpfile.txt)
[ "$OUT" = "test" ] || { rm -f tmpfile.txt; echo riker/coreutils 1; exit 1; }

# ls (just check that it runs and outputs something)
OUT=$("$COREUTILS_DIR/ls" .)
[ -n "$OUT" ] || { rm -f tmpfile.txt; echo riker/coreutils 1; exit 1; }

# Clean up
rm -f tmpfile.txt

echo coreutils/verify $?