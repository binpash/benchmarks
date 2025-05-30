#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
IN="$REPO_TOP/aurpkg/inputs/packages"
OUT="${OUT:-$REPO_TOP/aurpkg/outputs}"
KOALA_SHELL="${KOALA_SHELL:-bash}"
SCRIPT="./scripts/pacaur.sh"

for arg in "$@"; do
  if [ "$arg" = "--small" ]; then
    IN="$REPO_TOP/aurpkg/inputs/packages_small"
  elif [ "$arg" = "--min" ]; then
    IN="$REPO_TOP/aurpkg/inputs/packages_min"
  fi
done

mkdir -p "${OUT}"

# Set environment variables
BENCHMARK_CATEGORY="aurpkg"
export BENCHMARK_CATEGORY

BENCHMARK_SCRIPT="$(realpath "$SCRIPT")"
export BENCHMARK_SCRIPT

BENCHMARK_INPUT_FILE="$(realpath "$IN")"
export BENCHMARK_INPUT_FILE

echo "$SCRIPT"

if [ "$EUID" -eq 0 ]; then
  if ! id "user" &>/dev/null; then
    echo "Creating user 'user'..."
    useradd -m user
  fi

  echo "Running script as 'user'..."
  chown -R user:user "$OUT"
  $KOALA_SHELL "$BENCHMARK_SCRIPT" "$IN" "$OUT"

else
  echo "Not root, running script..."
  $KOALA_SHELL "$SCRIPT" "$IN" "$OUT"
fi

echo "$?"
