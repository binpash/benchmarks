#!/bin/bash
REPO_TOP=$(git rev-parse --show-toplevel)
eval_dir="${REPO_TOP}/vps-audit"
rm "${eval_dir}"/vps-audit-negate-report.txt
rm "${eval_dir}"/vps-audit-negate-processed.txt
rm "${eval_dir}"/vps-audit-report.txt
rm "${eval_dir}"/vps-audit-processed.txt
