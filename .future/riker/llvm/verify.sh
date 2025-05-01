#!/bin/sh
REPO_TOP="$(git rev-parse --show-toplevel)"
eval_dir="${REPO_TOP}/riker"
input_dir="${eval_dir}/input"
LLVM_BUILD_DIR="$input_dir/scripts/llvm/dev"
LLVM_BIN_DIR="$LLVM_BUILD_DIR/bin"
LLVM_AS="$LLVM_BIN_DIR/llvm-as"
LLVM_DIS="$LLVM_BIN_DIR/llvm-dis"
LLC="$LLVM_BIN_DIR/llc"

TMP_LL="test.ll"
TMP_BC="test.bc"
TMP_S="test.s"

cat <<EOF > "$TMP_LL"
; ModuleID = 'test'
source_filename = "test.c"
define i32 @main() {
entry:
  ret i32 0
}
EOF

"$LLVM_AS" "$TMP_LL" -o "$TMP_BC"
if [ $? -ne 0 ]; then
    echo riker/llvm 1
    exit 1
fi

"$LLVM_DIS" "$TMP_BC" -o /dev/null
if [ $? -ne 0 ]; then
    echo riker/llvm 1
    exit 1
fi

"$LLC" "$TMP_BC" -o "$TMP_S"
if [ $? -ne 0 ]; then
    echo riker/llvm 1
    exit 1
fi

"$LLVM_BIN_DIR/llvm-config" --version | grep -qE '^[0-9]+\.[0-9]+'
if [ $? -ne 0 ]; then
    echo riker/llvm 1
    exit 1
fi

rm -f "$TMP_LL" "$TMP_BC" "$TMP_S"

echo riker/llvm $?