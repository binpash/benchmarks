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
    pushd "${temp}" > /dev/null
    mkdir -p src/.git
    echo "echo This is a test" > src/startup.sh
}

tearDown() {
    popd > /dev/null
    rm -rf "${temp}"
}

testTarExtraOpts() {
    setupTests

    local tar_extra="--verbose --exclude .git"
    if "$SUT" --tar-extra "$tar_extra" src src.sh alabel startup.sh; then
        log_result "testTarExtraOpts" "PASS"
    else
        log_result "testTarExtraOpts" "FAIL" "Tar extra options failed."
        tearDown
        return 1
    fi

    tearDown
}

# Run the test
echo "Running testTarExtraOpts..."
if testTarExtraOpts; then
    log_result "testTarExtraOpts" "PASS"
else
    log_result "testTarExtraOpts" "FAIL"
fi

echo "Tests completed. Results logged in ${LOGFILE}"
