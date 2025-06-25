#!/bin/bash

cd "$(dirname "$0")" || exit 1

benchmarks="analytics bio ci-cd covid file-mod inference ml nlp oneliners pkg repl unixfun weather web-search web-search.new"

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
