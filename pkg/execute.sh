#!/bin/bash

TOP=$(realpath "$(dirname "$0")")


size=full
for arg in "$@"; do
    case "$arg" in
    --small) size=small ;;
    --min) size=min ;;
    esac
done

IN="$TOP/inputs/packages.$size"

KOALA_SHELL="${KOALA_SHELL:-bash}"
export BENCHMARK_CATEGORY="pkg"

SUITE_DIR="""$(realpath "$(dirname "$0")")"
export SUITE_DIR

export TIMEFORMAT=%R
cd "$SUITE_DIR" || exit 1

script_file="$TOP/scripts/pacaur.sh"

BENCHMARK_SCRIPT="$(realpath "$script_file")"
export BENCHMARK_SCRIPT

BENCHMARK_INPUT_FILE="$(realpath "$IN")"
export BENCHMARK_INPUT_FILE

echo "pacaur.sh"
OUT="$TOP/outputs/aurpkg.$size"
mkdir -p "${OUT}"
if [ "$EUID" -eq 0 ]; then
  if ! id "user" &>/dev/null; then
    echo "Creating user 'user'..."
    useradd -m user
  fi

  echo "Running script as 'user'..."
  chown -R user:user "$OUT"
  $KOALA_SHELL "$script_file" "$IN" "$OUT"

else
  echo "Not root, running script..."
  $KOALA_SHELL "$script_file" "$IN" "$OUT"
fi

echo "$?"

export INDEX="$TOP/inputs/index.$size.txt"

script_file="$TOP/scripts/proginf.sh"
export BENCHMARK_SCRIPT=$(realpath "$script_file")
export BENCHMARK_INPUT_FILE="$TOP/inputs/node_modules"

echo "proginf.sh"
$KOALA_SHELL "$script_file"
echo "$?"
