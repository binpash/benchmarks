#!/bin/bash

cd "$(dirname "$0")" || exit 1

size=min
for arg in "$@"; do
    case "$arg" in
        --small) size="--small" ;;
        --min)   size="--min" ;;
        --full)   size="" ;;
    esac
done

benchmarks="oneliners bio vps-audit file-enc log-analysis web-index media-conv covid-mts sklearn nlp aurpkg unix50 makeself riker max-temp"

log() { echo -e "[*] $1"; }

log "Running kick tires script"

log "Checking if user can run sudo"
if [[ $EUID -ne 0 ]] && ! sudo -v &> /dev/null; then
  echo "You need to have sudo privileges or be the root user to run this script (package installations will be required)."
  exit 1
fi

log "Running ./main.sh --min for each benchmark"

# Run ./main.sh --min in each benchmark directory
for bench in $benchmarks; do
  log "Running benchmark: $bench"
  ./main.sh "$size" "$bench" || {
    echo "Failed to run --min for $bench"
    exit 1
  }
done

log "-----------------------------------------------------"
log "Everything seems to be working fine. You can proceed!"
log "-----------------------------------------------------"
