#!/bin/bash

do_write_only() {
    # read one megabyte, write a whole gigabyte

    base="$(dirname "$0")/../1M.ones"
    ones="$(cat "$base")"

    i=0
    while [ "$i" -lt 1000 ]; do
        i=$((i + 1))
        printf "$ones" >> /dev/null
    done

}

do_write_only

sleep 1
