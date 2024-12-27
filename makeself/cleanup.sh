#! /bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
EVAL_DIR="${REPO_TOP}/makeself"
TEST_DIR="${EVAL_DIR}/test"
SHUNIT2_DIR="${TEST_DIR}/shunit2"

rm -rf "${SHUNIT2_DIR}"