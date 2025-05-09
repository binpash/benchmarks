#!/bin/bash

ASSIGNMENT_STRING=${1?No assignment was given to execute}
PRE_ENV_FILE=${2?No env file to run with given}
POST_EXEC_ENV=${3?No Riker env file given}

source "$PASH_TOP/compiler/orchestrator_runtime/speculative/pash_spec_init_setup.sh"

RUN=$(printf 'source %s; %s\n source ${RUNTIME_DIR}/pash_declare_vars.sh %s' "${PRE_ENV_FILE}" "${ASSIGNMENT_STRING}" "${POST_EXEC_ENV}")

bash -c "$RUN"
