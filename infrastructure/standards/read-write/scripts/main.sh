#!/bin/bash

# tests that we're counting read_chars and write_chars in child processes

do_copy() {
    dd if=/dev/zero of=/dev/null bs=1M count=250
}

do_copy

"$(dirname "$0")/sub.sh" 750

sleep 1
