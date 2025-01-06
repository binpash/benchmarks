#!/bin/bash

# IN="$1"
# OUT="$2"

# mkcd() { mkdir -p "$1" && cd "$1"; }

# # check if not running as root
# # test "$UID" -gt 0 || { info "don't run this as root!"; exit; }

# # set link to plaintext PKGBUILDs
# pkgbuild="https://aur.archlinux.org/cgit/aur.git/plain/PKGBUILD?h"

# run_tests() {
#     pgk=$1
#     mkcd "${OUT}/$pkg"

#     curl --insecure -o  PKGBUILD "$pkgbuild=$pkg" 2> /dev/null || echo ' '

#     #info "fetch required pgp keys from PKGBUILD"
#     #gpg --recv-keys $(sed -n "s:^validpgpkeys=('\([0-9A-Fa-fx]\+\)').*$:\1:p" PKGBUILD)
#     # Some failure is expected here, so we ignore the return code
#     makedeb -d >> ../$pkg.txt 2>&1
#     cd -
# }
# export -f run_tests

# # loop over required packages
# for pkg in $(cat ${IN} | tr '\n' ' ' ); 
# do  
#     echo "$pkg"
#     run_tests $pkg
# done

# Using GNU Parallel:

IN="$1"
OUT="$2"

mkcd() { mkdir -p "$1" && cd "$1"; }

# Set link to plaintext PKGBUILDs
pkgbuild="https://aur.archlinux.org/cgit/aur.git/plain/PKGBUILD?h"

run_tests() {
    pkg="$1"
    out_dir="$2"

    mkcd "${out_dir}/$pkg" || exit 1

    curl --insecure -o PKGBUILD "$pkgbuild=$pkg" 2> /dev/null || echo ' '

    # Info: Fetch required PGP keys from PKGBUILD (commented in the original script)
    # gpg --recv-keys $(sed -n "s:^validpgpkeys=('\([0-9A-Fa-fx]\+\)').*$:\1:p" PKGBUILD)
    makedeb -d >> "../$pkg.txt" 2>&1
    cd - > /dev/null || exit 1
}
export -f run_tests
export pkgbuild

# Read package names from the input file and process them in parallel
cat "$IN" | tr '\n' ' ' | parallel --jobs "$(nproc)" run_tests {} "$OUT"
