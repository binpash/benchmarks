#!/bin/bash
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="$REPO_TOP/tuft-weather"

for arg in "$@"; do
    case "$arg" in
        -f | --force)
            rm -rf $eval_dir/inputs/
            rm -rf $eval_dir/venv
            ;;
    esac
done
rm -rf $eval_dir/outputs/
rm -rf $eval_dir/plots/
rm -rf $eval_dir/full/
rm -rf $eval_dir/small/
rm -rf $eval_dir/min/