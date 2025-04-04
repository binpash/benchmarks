#!/bin/bash

REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/git-workflow"
outputs_dir="${eval_dir}/outputs"

if [ ! -d "$outputs_dir" ]; then
    echo "Outputs directory not found: $outputs_dir"
    exit 1
fi

for i in {1..5}; do
    status_log="${outputs_dir}/status-${i}.log"
    add_log="${outputs_dir}/add-${i}.log"
    commit_log="${outputs_dir}/commit-${i}.log"
    patchprep_log="${outputs_dir}/patchprep-${i}.log"
    
    for logfile in "$status_log" "$add_log" "$commit_log" "$patchprep_log"; do
        if [ ! -f "$logfile" ]; then
            echo "Missing log file: $logfile"
            exit 1
        elif [ ! -s "$logfile" ]; then
            echo "Log file is empty: $logfile"
            exit 1
        fi
    done
done
