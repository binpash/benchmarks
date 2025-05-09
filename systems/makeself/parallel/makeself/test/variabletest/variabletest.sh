#!/usr/bin/env bash
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

setupTests() {
    temp=$(mktemp -d -t XXXXX)
    cd "$temp"
    mkdir archive
    touch archive/file

    "$SUT" archive makeself-test.run "Test $1" echo "\\\"\${${1}}\\\""
}

testArchiveDir() {
    setupTests ARCHIVE_DIR
    local ans=$'./complicated\n dir\twith  spaces'
    mkdir "${ans}"
    mv ./makeself-test.run "${ans}/"

    local actual_archive_dir
    actual_archive_dir="$("${ans}/makeself-test.run" --quiet)"

    if [[ "${actual_archive_dir}" == "${ans}" ]]; then
        log_result "testArchiveDir" "PASS"
    else
        log_result "testArchiveDir" "FAIL" "Expected: ${ans}, Got: ${actual_archive_dir}"
    fi
    rm -rf "${temp}"
}

testTmpRoot() {
    setupTests TMPROOT
    local ans="${temp}"$'/complicated\n dir\twith  spaces'
    mkdir -p "${ans}"

    local actual_tmp_root
    actual_tmp_root="$(TMPDIR="${ans}" "./makeself-test.run" --quiet)"

    if [[ "${actual_tmp_root}" == "${ans}" ]]; then
        log_result "testTmpRoot" "PASS"
    else
        log_result "testTmpRoot" "FAIL" "Expected: ${ans}, Got: ${actual_tmp_root}"
    fi
    rm -rf "${temp}"
}

testUserPWD() {
    setupTests USER_PWD
    local ans="${temp}"$'/complicated\n dir\twith  spaces'
    mkdir -p "${ans}"
    cd "${ans}"

    local actual_user_pwd
    actual_user_pwd="$("${temp}/makeself-test.run" --quiet)"

    if [[ "${actual_user_pwd}" == "${ans}" ]]; then
        log_result "testUserPWD" "PASS"
    else
        log_result "testUserPWD" "FAIL" "Expected: ${ans}, Got: ${actual_user_pwd}"
    fi
    rm -rf "${temp}"
}

# Run tests
testArchiveDir
testTmpRoot
testUserPWD

echo "Tests completed. Results logged in ${LOGFILE}"
