#!/bin/bash

IN="$1"
OUT="$2"

mkcd() {
    mkdir -p "$1" || return 1
    cd "$1" || return 1
}

# set link to plaintext PKGBUILDs
URL='https://atlas.cs.brown.edu/data/aurpkg'

run_tests() {
    pkg=$1
    ORIG_DIR=$(pwd)

    mkcd "${OUT}/$pkg" || exit 1

    curl --insecure -o PKGBUILD "$URL/$pkg/PKGBUILD" 2>/dev/null || echo ' '

    #info "fetch required pgp keys from PKGBUILD"
    gpg --recv-keys $(sed -n "s:^validpgpkeys=('\([0-9A-Fa-fx]\+\)').*$:\1:p" PKGBUILD)
    # Some failure is expected here, so we ignore the return code
    makedeb -d >>"../$pkg.txt" 2>&1 || true
    cd "$ORIG_DIR" || exit 1
}

export -f run_tests

pkg_count=0

# loop over required packages
for pkg in $(cat ${IN}  | tr '\n' ' '); do
    pkg_count=$((pkg_count + 1))
    echo "$pkg"
    run_tests $pkg>"${OUT}/$pkg_count-$pkg.txt"
done
