#!/bin/bash

IN="$1"
OUT="$2"

mkcd() { mkdir -p "$1" && cd "$1"; }

# check if not running as root
# test "$UID" -gt 0 || { info "don't run this as root!"; exit; }

# set link to plaintext PKGBUILDs
pkgbuild="https://aur.archlinux.org/cgit/aur.git/plain/PKGBUILD?h"

run_tests() {
    pgk=$1
    mkcd "${OUT}/$pkg"

    curl --insecure -o  PKGBUILD "$pkgbuild=$pkg" 2> /dev/null || echo ' '

    #info "fetch required pgp keys from PKGBUILD"
    #gpg --recv-keys $(sed -n "s:^validpgpkeys=('\([0-9A-Fa-fx]\+\)').*$:\1:p" PKGBUILD)
    # Some failure is expected here, so we ignore the return code
    makedeb -d >> ../$pkg.txt 2>&1
    cd -
}
export -f run_tests

pkg_count=0

# loop over required packages
for pkg in $(cat ${IN} | tr '\n' ' ' ); 
do  
    pkg_count=$((pkg_count+1))
    echo "$pkg"
    run_tests $pkg>"${OUT}/$pkg_count.txt"
done
