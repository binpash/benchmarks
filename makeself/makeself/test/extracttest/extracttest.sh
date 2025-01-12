#!/usr/bin/env bash
set -eu
THIS="$(readlink -f "$0")"
THISDIR="$(dirname "${THIS}")"
SRCDIR="$(dirname "$(dirname "${THISDIR}")")"
SUT="${SRCDIR}/makeself.sh"
LOGFILE="${THISDIR}/test_results.log"

# Initialize the log file
echo "Test results:" > "${LOGFILE}"

log_result() {
    local test_name="$1"
    local result="$2"
    local details="${3:-}"
    echo "${result}: ${test_name} ${details}" >> "${LOGFILE}"
}

setupTests() {
    temp=$(mktemp -d -t XXXXX)
    cd "$temp"
    mkdir archive
    cp -a "$SRCDIR" archive/
    "$SUT" "$@" archive makeself-test.run "Test $*" echo Testing --tar-extra="--exclude .git"
}

run_test() {
    local test_name="$1"
    shift
    echo "Running ${test_name}..."
    if "$@"; then
        log_result "${test_name}" "PASS"
    else
        log_result "${test_name}" "FAIL"
    fi
}

testQuiet() {
    setupTests
    ./makeself-test.run --quiet
    return $?
}

testGzip() {
    setupTests --gzip
    ./makeself-test.run --check
    return $?
}

testBzip2() {
    setupTests --bzip2
    ./makeself-test.run --check
    return $?
}

testPBzip2() {
    if ! command -v pbzip2 >/dev/null 2>&1; then
        log_result "testPBzip2" "SKIPPED" "pbzip2 not available"
        return 0
    fi
    setupTests --pbzip2
    ./makeself-test.run --check
    return $?
}

testZstd() {
    if ! command -v zstd >/dev/null 2>&1; then
        log_result "testZstd" "SKIPPED" "zstd not available"
        return 0
    fi
    setupTests --zstd
    ./makeself-test.run --check
    return $?
}

# Run tests and log results
run_test "testQuiet" testQuiet
run_test "testGzip" testGzip
run_test "testBzip2" testBzip2
run_test "testPBzip2" testPBzip2
run_test "testZstd" testZstd

echo "Tests completed. Results logged in ${LOGFILE}"
