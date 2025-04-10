#!/bin/bash
set -e

REPO_TOP="$(git rev-parse --show-toplevel)"
PROTOBUF_DIR="$REPO_TOP/riker/input/scripts/protobuf/dev"
PROTOC_BIN="$PROTOBUF_DIR/src/protoc"
PROTOBUF_LIB="$PROTOBUF_DIR/src/.libs/libprotobuf.so"
INCLUDE_DIR="$PROTOBUF_DIR/src"

if [ ! -x "$PROTOC_BIN" ] || [ ! -f "$PROTOBUF_LIB" ]; then
  echo riker/protobuf 1
  exit 1
fi

TMP_DIR=$(mktemp -d) || { echo riker/protobuf 1; exit 1; }
cd "$TMP_DIR" || { echo riker/protobuf 1; exit 1; }

# Write test.proto
cat <<EOF > test.proto
syntax = "proto3";

message TestMessage {
  string test_field = 1;
}
EOF

# Compile proto to C++
"$PROTOC_BIN" --cpp_out=. test.proto > /dev/null 2>&1

if [ ! -f "test.pb.cc" ] || [ ! -f "test.pb.h" ]; then
  echo riker/protobuf 1
  cd .. >/dev/null 2>&1
  rm -rf "$TMP_DIR" >/dev/null 2>&1
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
  -lprotobuf > /dev/null 2>&1

# Run and check output
OUTPUT=$(./test_program 2>/dev/null)
EXPECTED_OUTPUT="Hello, Protobuf!"

if [ "$OUTPUT" != "$EXPECTED_OUTPUT" ]; then
  echo riker/protobuf 1
  cd .. >/dev/null 2>&1
  rm -rf "$TMP_DIR" >/dev/null 2>&1
  exit 1
fi

cd .. >/dev/null 2>&1
rm -rf "$TMP_DIR" >/dev/null 2>&1

echo riker/protobuf 0
exit 0
