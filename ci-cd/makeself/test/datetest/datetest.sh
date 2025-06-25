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

setUp() {
    temp="$(mktemp -dt datetest.XXXXX)"
    cd "${temp}"
    mkdir src
    echo "echo This is a test" > src/startup.sh
}

tearDown() {
    cd - > /dev/null
    rm -rf "${temp}"
}

testCurrentDate() {
    setUp
    "$SUT" src src.sh alabel startup.sh

    actual=$(strings src.sh | grep packaging || true)
    expected=$(LC_ALL=C date +"%b")

    if [[ ${actual} == *${expected}* ]]; then
        log_result "testCurrentDate" "PASS"
    else
        log_result "testCurrentDate" "FAIL" "Expected '${expected}', Found '${actual}'"
    fi
    tearDown
}

testDateSet() {
    setUp
    expected='Sat Mar  5 19:35:21 EST 2016'

    "$SUT" --packaging-date "${expected}" src src.sh alabel startup.sh

    actual=$(strings src.sh | grep "Date of packaging" || true)

    if [[ ${actual} == *${expected}* ]]; then
        log_result "testDateSet" "PASS"
    else
        log_result "testDateSet" "FAIL" "Expected '${expected}', Found '${actual}'"
    fi
    tearDown
}

testPackagingDateNeedsParameter() {
    setUp
    if ! "$SUT" --packaging-date src src.sh alabel startup.sh; then
        log_result "testPackagingDateNeedsParameter" "PASS"
    else
        log_result "testPackagingDateNeedsParameter" "FAIL" "Expected failure, but succeeded"
    fi
    tearDown
}

testByteforbyte() {
    setUp
    date='Sat Mar  3 19:35:21 EST 2016'

    "$SUT" --packaging-date "${date}" --tar-extra "--mtime 20160303" \
        src src.sh alabel startup.sh
    mv src.sh first
    "$SUT" --packaging-date "${date}" --tar-extra "--mtime 20160303" \
        src src.sh alabel startup.sh
    mv src.sh second

    if cmp first second > /dev/null; then
        log_result "testByteforbyte" "PASS"
    else
        log_result "testByteforbyte" "FAIL" "Files 'first' and 'second' differ"
    fi
    tearDown
}

# Run tests
testCurrentDate
testDateSet
testPackagingDateNeedsParameter
testByteforbyte

echo "Tests completed. Results logged in ${LOGFILE}"
