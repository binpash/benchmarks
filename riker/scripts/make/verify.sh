#!/bin/sh
set -e

# Path to the built GNU Make binary
REPO_TOP="$(git rev-parse --show-toplevel)"
MAKE_BIN="${REPO_TOP}/riker/input/scripts/make/dev/make"

# Verify binary exists
if [ ! -x "$MAKE_BIN" ]; then
  echo "make binary not found at $MAKE_BIN"
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
"$MAKE_BIN"

# Check output
if [ ! -f output.txt ]; then
  echo "make verification failed: output.txt not created"
  exit 1
fi

RESULT=$(cat output.txt | tr -d '[:space:]')
if [ "$RESULT" != "hellofrommake" ]; then
  echo "make verification failed: unexpected output '$RESULT'"
  exit 1
fi

# Clean up
cd ..
rm -rf "$WORKDIR"

echo make/verify $?
