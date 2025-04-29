#!/bin/bash
set -eu
THIS="$(readlink -f "$0")"
THISDIR="$(dirname "${THIS}")"
SUT="$(dirname "$(dirname "${THISDIR}")")/makeself.sh"
LOGFILE="${THISDIR}/test_results.log"

echo "Test results:" > "${LOGFILE}"

log_result() {
    local test_name="$1"
    local result="$2"
    local details="${3:-}"
    echo "${result}: ${test_name} ${details}" >> "${LOGFILE}"
}

withlsm() {
    local options="$@"
    (
        cd "${THISDIR}"
        mkdir -p lsmtest
        if ! "${SUT}" $options ./lsmtest ./lsmtest.run lsmtest ls -lah > /dev/null 2>&1; then
            log_result "withlsm" "FAIL" "Failed to create archive"
            return 1
        fi

        if ! ./lsmtest.run --lsm > "${THISDIR}/lsm_output.log" 2>&1; then
            log_result "withlsm" "FAIL" "Failed to extract archive or print LSM"
            rm -rf lsmtest lsmtest.run
            return 1
        fi

        log_result "withlsm" "PASS"
        rm -rf lsmtest lsmtest.run
    )
}

test_lsm_empty() {
    printf '' > lsm_empty.txt
    if withlsm --lsm lsm_empty.txt; then
        log_result "test_lsm_empty" "PASS"
    else
        log_result "test_lsm_empty" "FAIL"
    fi
    rm -f lsm_empty.txt
}

test_lsm_one_line() {
    printf 'one line\n' > lsm_one_line.txt
    if withlsm --lsm lsm_one_line.txt; then
        log_result "test_lsm_one_line" "PASS"
    else
        log_result "test_lsm_one_line" "FAIL"
    fi
    rm -f lsm_one_line.txt
}

test_lsm_one_line_without_nl() {
    printf 'one line without nl' > lsm_one_line_without_nl.txt
    if withlsm --lsm lsm_one_line_without_nl.txt; then
        log_result "test_lsm_one_line_without_nl" "PASS"
    else
        log_result "test_lsm_one_line_without_nl" "FAIL"
    fi
    rm -f lsm_one_line_without_nl.txt
}

# Run tests
echo "Running tests..."
test_lsm_empty
test_lsm_one_line
test_lsm_one_line_without_nl

echo "Tests completed. Results logged in ${LOGFILE}"
