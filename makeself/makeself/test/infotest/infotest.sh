#!/usr/bin/env bash
set -eu
THIS="$(readlink -f "$0")"
THISDIR="$(dirname "${THIS}")"
SRCDIR="$(dirname "$(dirname "${THISDIR}")")"
VERSION="$(cat "${SRCDIR}/VERSION")"
LOGFILE="${THISDIR}/test_results.log"

is_alpine_distro=false && [[ -f "/etc/alpine-release" ]] && is_alpine_distro=true

uncompressed_size="12 KB" && [[ $is_alpine_distro == true ]] && uncompressed_size="4 KB"

echo "Test results:" > "${LOGFILE}"

log_result() {
    local test_name="$1"
    local result="$2"
    local details="${3:-}"
    echo "${result}: ${test_name} ${details}" >> "${LOGFILE}"
}

# Take makeself options, generate a predefined archive, print --info to stdout.
#
# $@ : makeself options
haveInfo() {
    (
        cd "${SRCDIR}" || return 1
        mkdir -p infotest
        ./makeself.sh "$@" ./infotest ./infotest.run infotest ls -lah >/dev/null 2>&1
        if [[ $? -ne 0 ]]; then
            return 1
        fi
        ./infotest.run --info
        local rc=$?
        rm -rf infotest infotest.run
        return "${rc}"
    )
}

# Read want.info from stdin. Generate have.info using given options. Invoke
# diff want.info have.info and return its exit status
#
# $@ : makeself options
diffInfo() {
    local rc=""
    local valid_sizes="4 KB|12 KB" # Pipe-separated list of valid sizes for POSIX compliance
    cd "$(mktemp -d)" || return 1
    cat >want.info
    haveInfo "$@" >have.info

    # Replace the uncompressed size in have.info and want.info with a placeholder
    sed 's/Uncompressed size: .*$/Uncompressed size: VALID_SIZE/' have.info >have.modified.info
    sed 's/Uncompressed size: .*$/Uncompressed size: VALID_SIZE/' want.info >want.modified.info

    # Compare the modified files
    if diff want.modified.info have.modified.info >&2; then
        rc=0
    else
        rc=1
    fi

    # Extract the actual size from have.info
    actual_size=$(grep "Uncompressed size" have.info | awk '{print $3, $4}')

    # Check if the actual size is in the list of valid sizes
    if ! echo "$actual_size" | grep -Eq "^($valid_sizes)$"; then
        echo "Error: Uncompressed size '${actual_size}' is not valid. Expected one of: 4 KB or 12 KB" >&2
        rc=1
    fi

    rm -f have.info want.info have.modified.info want.modified.info
    return "${rc}"
}


# Run a test and log results
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

testDefault() {
    diffInfo --packaging-date "@0" <<EOF
Identification: infotest
Target directory: infotest
Uncompressed size: $uncompressed_size
Compression: gzip
Encryption: n
Date of packaging: @0
Built with Makeself version ${VERSION}
Build command was: ./makeself.sh \\
    "--packaging-date" \\
    "@0" \\
    "./infotest" \\
    "./infotest.run" \\
    "infotest" \\
    "ls" \\
    "-lah"
Script run after extraction:
     ls -lah
infotest will be removed after extraction
EOF
    return $?
}

testNocomp() {
    diffInfo --packaging-date "@0" --nocomp <<EOF
Identification: infotest
Target directory: infotest
Uncompressed size: $uncompressed_size
Compression: none
Encryption: n
Date of packaging: @0
Built with Makeself version ${VERSION}
Build command was: ./makeself.sh \\
    "--packaging-date" \\
    "@0" \\
    "--nocomp" \\
    "./infotest" \\
    "./infotest.run" \\
    "infotest" \\
    "ls" \\
    "-lah"
Script run after extraction:
     ls -lah
infotest will be removed after extraction
EOF
    return $?
}

testNotemp() {
    diffInfo --packaging-date "@0" --notemp <<EOF
Identification: infotest
Target directory: infotest
Uncompressed size: $uncompressed_size
Compression: gzip
Encryption: n
Date of packaging: @0
Built with Makeself version ${VERSION}
Build command was: ./makeself.sh \\
    "--packaging-date" \\
    "@0" \\
    "--notemp" \\
    "./infotest" \\
    "./infotest.run" \\
    "infotest" \\
    "ls" \\
    "-lah"
Script run after extraction:
     ls -lah
directory infotest is permanent
EOF
    return $?
}

# Run all tests and log results
run_test "testDefault" testDefault
run_test "testNocomp" testNocomp
run_test "testNotemp" testNotemp

echo "Tests completed. Results logged in ${LOGFILE}"
