#!/bin/sh
set -e
REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input"
AUTOCONF_BUILD_DIR="$input_dir/scripts/autoconf/dev"
AUTOCONF_BIN="$AUTOCONF_BUILD_DIR/bin/autoconf"
export PATH="$AUTOCONF_BUILD_DIR/bin:$PATH"
AUTOM4TE_PERLLIBDIR="$AUTOCONF_BUILD_DIR/share/autoconf"
export AUTOM4TE_PERLLIBDIR
status=0
if [ ! -x "$AUTOCONF_BIN" ]; then
    echo "Error: autoconf binary not found or not executable at $AUTOCONF_BIN"
    status=1
    echo riker/autoconf $status
    exit 1    
fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

cat > "$tmpdir/configure.ac" <<EOF
AC_INIT([testpackage], [1.0])
AC_CONFIG_SRCDIR([configure.ac])
AC_OUTPUT
EOF

cd "$tmpdir"
"$AUTOCONF_BIN"

if [ ! -x "./configure" ]; then
    status=1
    echo "Error: configure script not found or not executable."
    echo riker/autoconf $status
    exit 1
fi

./configure --help > /dev/null
echo riker/autoconf $?