#!/bin/bash

##
## This is the entry of pash-speculate and calls parallel-orch/orch.py
##

## Find the source code top directory
export PASH_SPEC_TOP=${PASH_SPEC_TOP:-$(realpath $(dirname $0))}
export PASH_TOP=${PASH_TOP:-$PASH_SPEC_TOP/deps/pash}

if [ -w /sys/fs/cgroup/ ]; then
    mkdir -p /sys/fs/cgroup/frontier
    total_mem=$(free | awk '/Mem:/ { print $2 }')
    protected_mem=$(python3 -c "print(int(${total_mem}*0.75) << 10)")
    chmod 666 /sys/fs/cgroup/cgroup.procs
    chmod 666 /sys/fs/cgroup/frontier/cgroup.procs
    if [ $(whoami) == "root" ]; then
	bash -c "echo $protected_mem > /sys/fs/cgroup/frontier/memory.min"
    fi
fi

if [ -n "$PASH_TMP_DIR" ]; then
    mkdir -p $PASH_TMP_DIR/tmp/pash_spec
    # echo $PASH_TMP_DIR >&2
    export PASH_SPEC_TMP_PREFIX="$(mktemp -d "$PASH_TMP_DIR/tmp/pash_spec/pash_XXXXXXX")"
else
    mkdir -p /tmp/pash_spec
    export PASH_SPEC_TMP_PREFIX="$(mktemp -d /tmp/pash_spec/pash_XXXXXXX)"
fi

# echo "PASH_SPEC_TMP_PREFIX: $PASH_SPEC_TMP_PREFIX" >&2

## Initialize the scheduler-server
export PASH_SPEC_SCHEDULER_SOCKET="${PASH_SPEC_TMP_PREFIX}/scheduler_socket"

## TODO: Replace this with a call to pa.sh (which will start the scheduler on its own).
# python3 "$PASH_SPEC_TOP/parallel-orch/orch.py" "$@"
"$PASH_TOP/pa.sh" --speculative "$@"
EXITCODE=$?

if [ -w /sys/fs/cgroup/ ]; then
    rmdir /sys/fs/cgroup/frontier
fi
exit $EXITCODE
