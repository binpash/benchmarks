#!/bin/bash

IN="$1"
OUT="$2"

mkcd() {
    mkdir -p "$1" || return 1
    cd "$1" || return 1
}

URL='https://atlas.cs.brown.edu/data/aurpkg'

grab_scalar() {
    awk -F= -v key="$1" '$1==key {sub(/^"|"$/, "", $2); print $2; exit}' PKGBUILD
}

grab_array() {
    awk -F'[()]' -v key="$1" '$1~key {sub(/^[[:space:]]+|[[:space:]]+$/, ""); print $2; exit}' PKGBUILD
}

run_tests() {
    pkg=$1
    ORIG_DIR=$(pwd)

    mkcd "${OUT}/$pkg" || exit 1

    curl --insecure -o PKGBUILD "$URL/$pkg/PKGBUILD" 2>/dev/null || echo ' '

    gpg --recv-keys $(sed -n "s:^validpgpkeys=('\([0-9A-Fa-fx]\+\)').*$:\1:p" PKGBUILD)
    {
        printf '[convert] %s -> pacscript\n' "$pkg"

        pkgname=$(grab_scalar pkgname)
        pkgver=$(grab_scalar pkgver)
        pkgdesc=$(grab_scalar pkgdesc)
        homepage=$(grab_scalar url)
        source=$(grab_array source)
        depends=$(grab_array depends)
        makedep=$(grab_array makedepends)
        sha256=$(awk '/^sha256sums=\(/{gsub(/'\''|\)/,""); print $2; exit}' PKGBUILD)

        [ -n "$pkgname" ] || pkgname=$pkg
        [ -n "$pkgver"  ] || pkgver=0.0.0
        [ -n "$pkgdesc" ] || pkgdesc='No description'
        [ -n "$homepage" ] || homepage="https://aur.archlinux.org/packages/$pkgname"

        cat >"${pkg}.pacscript" <<EOF
pkgname=$pkgname
pkgver=$pkgver
pkgdesc="$pkgdesc"
url=$homepage
source=($source)
sha256sums=('$sha256')
depends=($depends)
makedepends=($makedep)
arch=('amd64')

build() {
    cd "\$srcdir/\$pkgname-\$pkgver" || return 1
    make
}

package() {
    make DESTDIR="\$pkgdir" install
}
EOF
        printf '[convert] pacscript saved as %s.pacscript\n' "$pkg"
    } > "../${pkg}.txt" 2>&1 || echo "Failed to convert $pkg"

    cd "$ORIG_DIR" || exit 1
}

export -f run_tests
export OUT

parallel --line-buffer 'echo {}; run_tests {}' ::: $(cat "$IN")
