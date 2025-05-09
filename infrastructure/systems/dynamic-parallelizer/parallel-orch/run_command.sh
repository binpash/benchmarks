#!/bin/bash

## TODO: Not sure if it is OK and ideal to give this a string
export CMD_STRING=${1?No command was given to execute}
export TRACE_FILE=${2?No trace file path given}
export STDOUT_FILE=${3?No stdout file given}
export LATEST_ENV_FILE=${4?No env file to run with given}
export SANDBOX_DIR=${5?No sandbox dir given}
export TMPDIR=${6?No tmp dir given}
export EXEC_MODE=${7?No execution mode given}
export CMD_ID=${8?No command id given}
export POST_EXEC_ENV=${9?No Riker env file given}
export EXECUTION_ID=${10?No execution id given}
LOWER_DIRS=${11?No lower dirs}

## KK 2023-04-24: Not sure this should be run every time we run a command
## GL 2023-07-08: Tests seem to pass without it
source "$PASH_TOP/compiler/orchestrator_runtime/speculative/pash_spec_init_setup.sh"

if [ "standard" == "$EXEC_MODE" ]; then
    echo $$ > /sys/fs/cgroup/frontier/cgroup.procs
elif [ "speculate" == "$EXEC_MODE" ]; then
    renice 20 -p $$ >/dev/null
fi

# mkdir -p /tmp/pash_spec/a
# mkdir -p /tmp/pash_spec/b
# export SANDBOX_DIR="$(mktemp -d /tmp/pash_spec/a/sandbox_XXXXXXX)/"
# export TEMPDIR="$(mktemp -d /tmp/pash_spec/b/sandbox_XXXXXXX)"
# echo tempdir $TEMPDIR
# echo sandbox $SANDBOX_DIR

${RUNTIME_LIBRARY_DIR}/fd_util -f "${LATEST_ENV_FILE}.fds" -p ${STDOUT_FILE} bash "${PASH_SPEC_TOP}/deps/try/try" -D "${SANDBOX_DIR}" -L "${LOWER_DIRS}" "${PASH_SPEC_TOP}/parallel-orch/template_script_to_execute.sh"
exit_code=$?
## Only used for debugging
# ls -R "${SANDBOX_DIR}/upperdir" 1>&2
# out=`head -3 $SANDBOX_DIR/upperdir/$TRACE_FILE`
## Send a message to the scheduler socket
## Assumes "${PASH_SPEC_SCHEDULER_SOCKET}" is set and exported

## Pass the proper exit code
msg="CommandExecComplete:${CMD_ID}|Exec id:${EXECUTION_ID}|Sandbox dir:${SANDBOX_DIR}|Trace file:${TRACE_FILE}|Tempdir:${TEMPDIR}"
daemon_response=$(pash_spec_communicate_scheduler_just_send "$msg") # Blocking step, daemon will not send response until it's safe to continue
(exit $exit_code)
