#!/usr/bin/env bash
set -eu
THIS="$(readlink -f "$0")"
THISDIR="$(dirname "${THIS}")"
SRCDIR="$(dirname "$(dirname "${THISDIR}")")"
SUT="${SRCDIR}/makeself.sh"
LOGFILE="${THISDIR}/test_results.log"

echo "Test results:" > "${LOGFILE}"

log_result() {
    local test_name="$1"
    local result="$2"
    local details="${3:-}"
    echo "${result}: ${test_name} ${details}" >> "${LOGFILE}"
}

setupTests() {
    local temp_path
    temp_path="$(mktemp -dt appendtest.XXXXXX)"
    cd "$temp_path"
    mkdir -p archive
    cp -a "$SRCDIR" archive/
    if ! "$SUT" "$@" archive makeself-test.run "Test $*" echo Testing --tar-extra="--exclude .git"; then
        log_result "setupTests" "FAIL" "Failed to create archive"
        exit 1
    fi
}

testExtraBytes() {
    setupTests --sha256

    if ./makeself-test.run --check; then
        log_result "testExtraBytes: Initial Check" "PASS"
    else
        log_result "testExtraBytes: Initial Check" "FAIL"
        return 1
    fi

    echo "Adding a bunch of random characters at the end!!" >> makeself-test.run

    if ! ./makeself-test.run --check; then
        log_result "testExtraBytes: Corrupted Archive Check" "PASS"
    else
        log_result "testExtraBytes: Corrupted Archive Check" "FAIL"
        return 1
    fi
}

testTruncated() {
    setupTests --sha256

    if ./makeself-test.run --check; then
        log_result "testTruncated: Initial Check" "PASS"
    else
        log_result "testTruncated: Initial Check" "FAIL"
        return 1
    fi

    dd if=makeself-test.run of=truncated.run bs=1 count=34303

    if ! truncated.run --check; then
        log_result "testTruncated: Truncated Archive Check" "PASS"
    else
        log_result "testTruncated: Truncated Archive Check" "FAIL"
        return 1
    fi
}

# Run tests
echo "Running testExtraBytes..."
if testExtraBytes; then
    log_result "testExtraBytes" "PASS"
else
    log_result "testExtraBytes" "FAIL"
fi

echo "Running testTruncated..."
if testTruncated; then
    log_result "testTruncated" "PASS"
else
    log_result "testTruncated" "FAIL"
fi

echo "Tests completed. Results logged in ${LOGFILE}"
