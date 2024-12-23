#!/bin/bash
set -e

REPO_TOP=$(git rev-parse --show-toplevel)
EVAL_DIR="${REPO_TOP}/makeself"
TEST_DIR="${EVAL_DIR}/test"
SHUNIT2_REPO="https://github.com/kward/shunit2.git"
SHUNIT2_DIR="${TEST_DIR}/shunit2"
SHUNIT2_COMMIT="47be8b23a46a7897e849f1841f0fb704d34d0f6e"

if [[ ! -d "$EVAL_DIR" ]]; then
    echo "Cloning makeself repository..."
    git clone https://github.com/megastep/makeself.git "$EVAL_DIR"
fi

cd "$TEST_DIR"

if [[ -d "$SHUNIT2_DIR" && -z "$(ls -A "$SHUNIT2_DIR")" ]]; then
    echo "Removing empty shunit2 directory..."
    rm -rf "$SHUNIT2_DIR"
    echo "Cloning shunit2 repository..."
    git clone "$SHUNIT2_REPO" "$SHUNIT2_DIR"
    cd "$SHUNIT2_DIR"
    git checkout "$SHUNIT2_COMMIT"
fi

cd "$EVAL_DIR"
# Run the test suite
echo "Running makeself test suite..."
if make test; then
    echo "All tests passed successfully."
else
    echo "Some tests failed. Check the output above for details."
    exit 1
fi
