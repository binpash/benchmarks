#!/bin/bash

mkdir -p $2
pure_func() {
    tempfile=$(mktemp)
    cat > $tempfile
    # extract DNS queries
    tcpdump -nn -r $tempfile -A 'port 53' 2> /dev/null | sort | uniq |grep -Ev '(com|net|org|gov|mil|arpa)' 2> /dev/null
    # extract URL
    tcpdump -nn -r $tempfile -s 0 -v -n -l 2> /dev/null | egrep -i "POST /|GET /|Host:" 2> /dev/null
    # extract passwords
    tcpdump -nn -r $tempfile -s 0 -A -n -l 2> /dev/null | egrep -i "POST /|pwd=|passwd=|password=|Host:" 2> /dev/null
    # extract telnet login/password
    tcpdump -nn -r $tempfile -s 0 -A -n -l 'port 23' 2> /dev/null | egrep -i "login:|password:" 2> /dev/null

    rm -f $tempfile
}
export -f pure_func

for item in $1/*;
do
    logname=$2/$(basename $item).log;
    cat $item | pure_func > $logname
done
