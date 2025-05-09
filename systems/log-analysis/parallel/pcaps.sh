#!/bin/bash

mkdir -p "$2"

pure_func() {
    tempfile=$(mktemp)
    cat > $tempfile 
    tcpdump -nn -r $tempfile -A 'port 53' 2> /dev/null | sort | uniq |grep -Ev '(com|net|org|gov|mil|arpa)' 2> /dev/null
    # extract URL
    tcpdump -nn -r $tempfile -s 0 -v -n -l 2> /dev/null | egrep -i "POST /|GET /|Host:" 2> /dev/null
    # extract passwords
    tcpdump -nn -r $tempfile -s 0 -A -n -l 2> /dev/null | egrep -i "POST /|pwd=|passwd=|password=|Host:" 2> /dev/null

    rm -f $tempfile
} # parallelizing the inner workload yielded worse performance

export -f pure_func

find "$1" -type f | parallel "cat {} | pure_func > $2/{/}.log"
