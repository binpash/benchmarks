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

testNoCheck() {
    local archive_dir file_name

    # Create a directory with a simple payload.
    archive_dir="$(mktemp -dt archive_dir.XXXXXX)"
    (
        cd "${archive_dir}"
        touch foo.txt bar.txt qux.txt
    )

    # Create a self-extracting archive.
    file_name="$(mktemp -t file_name.XXXXXX)"
    if "${SUT}" --nox11 --sha256 "${archive_dir}" "${file_name}" "no check test" true; then
        log_result "testNoCheck: Create Archive" "PASS"
    else
        log_result "testNoCheck: Create Archive" "FAIL"
        rm -rf "${archive_dir}" "${file_name}"
        return 1
    fi

    # Archive verification enabled
    printf '\nArchive verification enabled:\n' >&2
    sync
    if "${file_name}" 2>&1 | grep -qF 'Verifying archive integrity...'; then
        log_result "testNoCheck: Verify Archive (Enabled)" "PASS"
    else
        log_result "testNoCheck: Verify Archive (Enabled)" "FAIL"
        rm -rf "${archive_dir}" "${file_name}"
        return 1
    fi

    # Archive verification disabled
    printf '\nArchive verification disabled:\n' >&2
    if SETUP_NOCHECK=1 "${file_name}" 2>&1 | grep -qFv 'Verifying archive integrity...'; then
        log_result "testNoCheck: Verify Archive (Disabled)" "PASS"
    else
        log_result "testNoCheck: Verify Archive (Disabled)" "FAIL"
        rm -rf "${archive_dir}" "${file_name}"
        return 1
    fi

    # Clean up
    rm -rf "${archive_dir}" "${file_name}"
}

# Run the test
echo "Running testNoCheck..."
if testNoCheck; then
    log_result "testNoCheck" "PASS"
else
    log_result "testNoCheck" "FAIL"
fi

echo "Tests completed. Results logged in ${LOGFILE}"
