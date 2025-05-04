#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input/scripts"
scripts_dir="${eval_dir}/scripts"

KOALA_SHELL=${KOALA_SHELL:-bash}
(cd "$input_dir/lsof/dev" && $KOALA_SHELL "$scripts_dir/lsof/build.sh")

