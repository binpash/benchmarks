#!/bin/sh
set -e
REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input"
CALC_BUILD_DIR="$input_dir/scripts/calc/dev"
CALC_BIN="$CALC_BUILD_DIR/bin"
CALC_BIN="./checkout/calc"

if [ ! -x "$CALC_BIN" ]; then
  echo riker/calc 1
  exit 1
fi

RAW_RESULT=$("$CALC_BIN" 3 + 4)
RESULT=$(echo "$RAW_RESULT" | tr -d '[:space:]')

if [ "$RESULT" != "7" ]; then
  echo riker/calc 1
  exit 1
fi

echo calc $?