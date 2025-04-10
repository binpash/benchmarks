#!/bin/bash
set -e

REPO_TOP="$(git rev-parse --show-toplevel)"
PROTOBUF_DIR="$REPO_TOP/riker/input/scripts/protobuf/dev"
PROTOC_BIN="$PROTOBUF_DIR/src/protoc"
PROTOBUF_LIB="$PROTOBUF_DIR/src/.libs/libprotobuf.so"
INCLUDE_DIR="$PROTOBUF_DIR/src"

# Check compiler and library
if [ ! -x "$PROTOC_BIN" ]; then
  echo "Protobuf verification failed: protoc compiler not found at $PROTOC_BIN"
  exit 1
fi

if [ ! -f "$PROTOBUF_LIB" ]; then
  echo "Protobuf verification failed: protobuf library not found at $PROTOBUF_LIB"
  exit 1
fi

# Create a temp build directory
TMP_DIR=$(mktemp -d)
cd "$TMP_DIR"

# Write test.proto
cat <<EOF > test.proto
syntax = "proto3";

message TestMessage {
  string test_field = 1;
}
EOF

# Compile proto to C++
"$PROTOC_BIN" --cpp_out=. test.proto

if [ ! -f "test.pb.cc" ] || [ ! -f "test.pb.h" ]; then
  echo "Protobuf verification failed: Generated source files are missing"
  exit 1
fi

# Write test.cpp
cat <<EOF > test.cpp
#include "test.pb.h"
#include <iostream>

int main() {
  TestMessage msg;
  msg.set_test_field("Hello, Protobuf!");
  std::cout << msg.test_field() << std::endl;
  return 0;
}
EOF

# Compile test program with proper include and lib path
g++ test.cpp test.pb.cc -o test_program \
  -I. -I"$INCLUDE_DIR" \
  -L"$(dirname "$PROTOBUF_LIB")" \
  -Wl,-rpath,"$(dirname "$PROTOBUF_LIB")" \
  -lprotobuf

# Run and check output
OUTPUT=$(./test_program)
EXPECTED_OUTPUT="Hello, Protobuf!"

if [ "$OUTPUT" != "$EXPECTED_OUTPUT" ]; then
  echo riker/protobuf 1
  cd ..
  rm -rf "$TMP_DIR"
  exit 1
fi

cd ..
rm -rf "$TMP_DIR"

echo riker/protobuf $?
