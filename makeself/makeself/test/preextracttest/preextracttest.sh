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

setUp() {
    temp=$(mktemp -d -t XXXXX)
    pushd "${temp}" > /dev/null
    mkdir src
    echo "echo This is a test" > src/startup.sh
    chmod a+x src/startup.sh
}

tearDown() {
    popd > /dev/null
    rm -rf "${temp}"
}

testPreextractOpts() {
    setUp
    echo 'echo A complex pre-extraction script.
    sleep 99 &
    cat a.txt 2>/dev/null || cat b.txt && cat c.txt
    echo "$$ Some\toutput\n\a\b\0777 $var1 ${var2} `cat var3.txt` $(env)" > text.txt
  ' > preextract.sh

    "$SUT" --nox11 --preextract preextract.sh src src.sh alabel ./startup.sh
    if [[ $? -eq 0 ]]; then
        log_result "testPreextractOpts: Create Archive" "PASS"
    else
        log_result "testPreextractOpts: Create Archive" "FAIL"
        tearDown
        return 1
    fi

    ./src.sh --show-preextract > show-preextract.out
    if diff preextract.sh show-preextract.out > /dev/null; then
        log_result "testPreextractOpts: Verify Preextract" "PASS"
    else
        log_result "testPreextractOpts: Verify Preextract" "FAIL"
        tearDown
        return 1
    fi

    tearDown
}

testWithNoPreextractOpts() {
    setUp
    "$SUT" src src.sh alabel ./startup.sh
    if ./src.sh --show-preextract; then
        log_result "testWithNoPreextractOpts" "FAIL"
    else
        log_result "testWithNoPreextractOpts" "PASS"
    fi
    tearDown
}

testPreextractRun() {
    setUp
    echo 'echo Validating provided options...' > preextract.sh
    "$SUT" --nox11 --preextract preextract.sh src src.sh alabel ./startup.sh
    if ./src.sh | grep -qF 'Validating provided options...'; then
        log_result "testPreextractRun" "PASS"
    else
        log_result "testPreextractRun" "FAIL"
    fi
    tearDown
}

testPreextractNoexec() {
    setUp
    echo 'exit 2' > preextract.sh
    "$SUT" --preextract preextract.sh src src.sh alabel ./startup.sh
    if ./src.sh --noexec; then
        log_result "testPreextractNoexec" "PASS"
    else
        log_result "testPreextractNoexec" "FAIL"
    fi
    tearDown
}

testPreextractArgs() {
    setUp
    echo 'echo $*' > preextract.sh
    "$SUT" --nox11 --preextract preextract.sh src src.sh alabel ./startup.sh --logdir /var/log
    test_cmd='./src.sh -- --env dev'

    if eval "${test_cmd}" | grep -qF -- '--logdir /var/log --env dev'; then
        log_result "testPreextractArgs" "PASS"
    else
        log_result "testPreextractArgs" "FAIL"
    fi
    tearDown
}

testPreextractEnvPassing() {
    setUp
    echo 'echo "export INSTALLATION_DIR=/usr/bin" > preextract.env' > preextract.sh
    echo '. ./preextract.env; echo $INSTALLATION_DIR' > src/startup.sh
    "$SUT" --nox11 --preextract preextract.sh src src.sh alabel ./startup.sh
    if ./src.sh | grep -qF '/usr/bin'; then
        log_result "testPreextractEnvPassing" "PASS"
    else
        log_result "testPreextractEnvPassing" "FAIL"
    fi
    tearDown
}

# Run all tests
testPreextractOpts
testWithNoPreextractOpts
testPreextractRun
testPreextractNoexec
testPreextractArgs
testPreextractEnvPassing

echo "Tests completed. Results logged in ${LOGFILE}"
