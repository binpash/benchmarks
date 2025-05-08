#!/bin/bash

benchmarks="oneliners bio vps-audit file-enc log-analysis web-index media-conv covid-mts sklearn nlp aurpkg unix50 makeself riker max-temp"

log() { echo -e "[*] $1"; }

log "Running kick tires script"

log "Checking if user can run sudo"
if [[ $EUID -ne 0 ]] && ! sudo -v &> /dev/null; then
  echo "You need to have sudo privileges or be the root user to run this script (package installations will be required)."
  exit 1
fi


log "Running ./main.sh --min for each benchmark"

# STEP 1: Run ./main.sh --min in each benchmark directory
for bench in $benchmarks; do
  log "Running benchmark: $bench"
  ./main.sh --min "$bench" || {
    echo "Failed to run --min for $bench"
    exit 1
  }
done
