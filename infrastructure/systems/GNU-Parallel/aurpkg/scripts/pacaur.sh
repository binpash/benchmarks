#!/bin/bash
# Using GNU Parallel:
IN="$1"
OUT="$2"

mkcd() { mkdir -p "$1" && cd "$1"; }

# check if not running as root
# test "$UID" -gt 0 || { info "don't run this as root!"; exit; }

# Set link to plaintext PKGBUILDs
pkgbuild="https://aur.archlinux.org/cgit/aur.git/plain/PKGBUILD?h"

run_tests() {
    pkg="$1"
    echo "$pkg"
    mkcd "${OUT}/$pkg"

    curl --insecure -o PKGBUILD "$pkgbuild=$pkg" 2>/dev/null || echo ' '

    # Fetch required pgp keys from PKGBUILD (optional)
    # gpg --recv-keys $(sed -n "s:^validpgpkeys=('\([0-9A-Fa-fx]\+\)').*$:\1:p" PKGBUILD)
    # Some failure is expected here, so we ignore the return code
    makedeb -d >> "../$pkg.txt" 2>&1
    cd - 
}
export -f run_tests mkcd
export pkgbuild

# Read package names from the input file and process them in parallel
parallel run_tests :::: "$IN"
