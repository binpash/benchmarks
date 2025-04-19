#!/usr/bin/env bash
set -eu

THIS="$(readlink -f "$0")"
THISDIR="$(dirname "${THIS}")"
SRCDIR="$(dirname "$(dirname "${THISDIR}")")"
SUT="${SRCDIR}/makeself.sh"
WHAT="$(basename "${THIS}")"
LOGFILE="${THISDIR}/test_results.log"
echo "Test results:" > "${LOGFILE}"

log_result() {
    local test_name="$1"
    local result="$2"
    local details="${3:-}"
    echo "${result}: ${test_name} ${details}" >> "${LOGFILE}"
}

WORK_DIR=$(mktemp -d -t makeself_test.XXXXXX)
trap 'rm -rf "${WORK_DIR}"' EXIT 

readonly archive_dir_create="${WORK_DIR}/archive_dir_create"
readonly archive_dir_append="${WORK_DIR}/archive_dir_append"
mkdir -p "${archive_dir_create}" "${archive_dir_append}"
touch "${archive_dir_create}/fee" "${archive_dir_create}/fie"
touch "${archive_dir_append}/foe" "${archive_dir_append}/fum"

evalAssert() {
    eval "$@"
    if [[ $? -eq 0 ]]; then
        return 0
    else
        return 1
    fi
}

# $1 : file_name
doInfoListCheckExec() {
    evalAssert "$1 --info" &&
    evalAssert "$1 --list" &&
    evalAssert "$1 --check" &&
    evalAssert "$1"
}

# $1 : file_name
# rest : content basenames
assertFileContains() {
    local file_name="$(readlink -f "$1")"
    shift
    local target="${file_name}.d"
    rm -rf "${target}"
    mkdir -p "${target}"
    evalAssert "${file_name} --target ${target}" || return 1

    local expected
    expected=$(echo -e "$@" | sort)
    local actual
    actual=$(find "${target}" -type f -exec basename -a {} + | sort)

    if [[ "${actual}" == "${expected}" ]]; then
        rm -rf "${target}"
        return 0
    else
        rm -rf "${target}"
        return 1
    fi
}

# $@ : makeself options
doTestOpts() {
    local stem
    stem="$(printf '%s' "${WHAT}" "$@" | tr -sc '[:alnum:]_.-' '_')"
    local file_name="${WORK_DIR}/${stem}.run"

    if ! evalAssert "${SUT}" "$@" --sha256 "${archive_dir_create}" "${file_name}" "${stem}" "echo ${stem}"; then
        log_result "${stem}" "FAIL" "Error creating archive"
        return 1
    fi

    file_name="$(readlink -f "${file_name}")"
    doInfoListCheckExec "${file_name}" || {
        log_result "${stem}" "FAIL" "Error in info/list/check execution"
        return 1
    }
    assertFileContains "${file_name}" "fee\nfie" || {
        log_result "${stem}" "FAIL" "Archive content mismatch after creation"
        return 1
    }

    if ! evalAssert "${SUT}" "$@" --sha256 --append "${archive_dir_append}" "${file_name}"; then
        log_result "${stem}" "FAIL" "Error appending archive"
        return 1
    fi

    doInfoListCheckExec "${file_name}" || {
        log_result "${stem}" "FAIL" "Error in info/list/check execution after append"
        return 1
    }
    assertFileContains "${file_name}" "fee\nfie\nfoe\nfum" || {
        log_result "${stem}" "FAIL" "Archive content mismatch after append"
        return 1
    }

    rm -f "${file_name}"
    log_result "${stem}" "PASS"
}

# $1 : compression option
doTestComp() {
    if ! command -v "${1#--*}" >/dev/null 2>&1; then
        log_result "doTestComp ${1}" "SKIPPED" "Missing command: ${1#--*}"
        return 0
    fi
    doTestOpts "$1"
}

testDefault() { doTestOpts; }

testNocomp() { doTestOpts --nocomp; }

testBase64() { doTestComp --base64; }
testBzip2() { doTestComp --bzip2; }
testCompress() { doTestComp --compress; }
testGzip() { doTestComp --gzip; }
testLz4() { doTestComp --lz4; }
testLzo() { doTestComp --lzo; }
testPbzip2() { doTestComp --pbzip2; }
testPigz() { doTestComp --pigz; }
testXz() { doTestComp --xz; }
testZstd() { doTestComp --zstd; }

# Run all tests
cd "${WORK_DIR}"
testDefault
testNocomp
testBase64
testBzip2
testCompress
testGzip
testLz4
testLzo
testPbzip2
testPigz
testXz
testZstd

echo "Tests completed. Results logged in ${LOGFILE}"
