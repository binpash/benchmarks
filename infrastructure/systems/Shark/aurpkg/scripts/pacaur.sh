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

# loop over required packages
# for pkg in $(cat ${IN} | tr '\n' ' ' ); 
# do  
#     echo "$pkg"
#     run_tests $pkg
# done

for pkg in $(tr '\n' ' ' < "$IN"); do
    echo "Processing $pkg" &
    (run_tests $pkg) &
done
wait

#!/bin/bash

# IN="$1"
# OUT="$2"

# mkcd() { mkdir -p "$1"; }

# # Fetch and test packages
# run_tests() {
#     pkg="$1"
#     pkg_dir="${OUT}/$pkg"
#     mkcd "$pkg_dir"

#     # Fetch PKGBUILD only if not already downloaded
#     [[ -f "${pkg_dir}/PKGBUILD" ]] || curl --insecure -o "${pkg_dir}/PKGBUILD" "${pkgbuild}=$pkg" 2>/dev/null || {
#         echo "Failed to fetch PKGBUILD for $pkg"
#         return
#     }

#     # Run makedeb and save output
#     makedeb -d > "${pkg_dir}/${pkg}.txt" 2>&1
# }
# export -f run_tests

# MAX_JOBS=$(nproc) # Adjust based on system capacity
# current_jobs=0

# for pkg in $(<"$IN"); do
#     run_tests "$pkg" &
#     ((current_jobs++))

#     # Wait if the number of jobs reaches the limit
#     if ((current_jobs >= MAX_JOBS)); then
#         wait -n  # Wait for any job to complete
#         ((current_jobs--))
#     fi
# done

# # Wait for any remaining jobs to finish
# wait
