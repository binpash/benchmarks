#!/bin/bash

cd "$(dirname "$0")" || exit 1

benchmarks="aurpkg bio covid-mts dpt file-enc git-workflow llm log-analysis makeself max-temp media-conv nlp oneliners prog-inf riker sklearn teraseq tuft-weather unix50 vps-audit web-index"

log() { echo -e "[*] $1"; }

log "Running kick tires script"

# log "Checking if user can run sudo"
# if [[ $EUID -ne 0 ]] && ! sudo -v &> /dev/null; then
#   echo "You need to have sudo privileges or be the root user to run this script (package installations will be required)."
#   exit 1
# fi

log "Running ./main.sh for each benchmark"

# Run ./main.sh for each benchmark
for bench in $benchmarks; do
  log "Running benchmark: $bench"
  ./main.sh "$@" "$bench" || {
    echo "Failed to run $bench"
    exit 1
  }
done

log "-----------------------------------------------------"
log "Everything seems to be working fine. You can proceed!"
log "-----------------------------------------------------"
