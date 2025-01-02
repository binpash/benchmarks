#!/bin/bash

REPO_TOP="$(git rev-parse --show-toplevel)"

benches=$(python3 "$REPO_TOP/infrastructure/all_scripts.py" | sort)

for bench in $benches; do
    bash $REPO_TOP/$bench/deps.sh
done

for bench in $benches; do
    bash $REPO_TOP/$bench/input.sh
done

for bench in $benches; do
    python3 $REPO_TOP/infrastructure/run_dynamic.py $bench
done

touch "$REPO_TOP/infrastructure/target/collect_dynamic_logs.touch"
