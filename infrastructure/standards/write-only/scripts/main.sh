#!/bin/bash

do_write_only() {
    # read one megabyte, write a whole gigabyte

    base="$(dirname "$0")/../1M.ones"
    ones="$(cat "$base")"

    for ((i=0; i < 1000; i++)) {
        printf "$ones" >> /dev/null
    }

}

do_write_only

sleep 1
