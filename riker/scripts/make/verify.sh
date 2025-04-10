#!/bin/sh
set -e

# Path to the built GNU Make binary
REPO_TOP="$(git rev-parse --show-toplevel)"
MAKE_BIN="${REPO_TOP}/riker/input/scripts/make/dev/make"

# Verify binary exists
if [ ! -x "$MAKE_BIN" ]; then
  echo riker/make 1
  exit 1
fi

# Create a temporary test project
WORKDIR="$(mktemp -d)"
cd "$WORKDIR"

# Create a minimal Makefile and target
cat <<'EOF' > Makefile
all:
	echo "hello from make" > output.txt
EOF

# Run make
"$MAKE_BIN" > /dev/null 2>&1

# Check output
if [ ! -f output.txt ]; then
  echo riker/make 1
  exit 1
fi

RESULT=$(cat output.txt | tr -d '[:space:]')
if [ "$RESULT" != "hellofrommake" ]; then
  echo riker/make 1
  cd ..
  rm -rf "$WORKDIR"
  exit 1
fi

# Clean up
cd ..
rm -rf "$WORKDIR"

echo riker/make 0
exit 0
