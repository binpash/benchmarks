#!/bin/bash

# IN=${IN:-/dependency_untangling/pcap_data}
# OUT=${OUT:-$PASH_TOP/evaluation/distr_benchmarks/dependency_untangling/input/output/pcap-analysis}
# mkdir -p $2
# pure_func() {
#     tempfile=$(mktemp)
#     cat > $tempfile 
#     tcpdump -nn -r $tempfile -A 'port 53' 2> /dev/null | sort | uniq |grep -Ev '(com|net|org|gov|mil|arpa)' 2> /dev/null
#     # extract URL
#     tcpdump -nn -r $tempfile -s 0 -v -n -l 2> /dev/null | egrep -i "POST /|GET /|Host:" 2> /dev/null
#     # extract passwords
#     tcpdump -nn -r $tempfile -s 0 -A -n -l 2> /dev/null | egrep -i "POST /|pwd=|passwd=|password=|Host:" 2> /dev/null

#     rm -f $tempfile
# }
# export -f pure_func

# for item in $1/*;
# do
#     logname=$2/$(basename $item).log;
#     cat $item | pure_func > $logname
# done

mkdir -p "$2"

pure_func() {
    tcpdump -nn -r /dev/stdin -A 'port 53' 2> /dev/null | sort | uniq | grep -Ev '(com|net|org|gov|mil|arpa)' 2> /dev/null
    # Extract URL
    tcpdump -nn -r /dev/stdin -s 0 -v -n -l 2> /dev/null | egrep -i "POST /|GET /|Host:" 2> /dev/null
    # Extract passwords
    tcpdump -nn -r /dev/stdin -s 0 -A -n -l 2> /dev/null | egrep -i "POST /|pwd=|passwd=|password=|Host:" 2> /dev/null
}
export -f pure_func

for item in "$1"/*; do
    logname="$2/$(basename "$item").log"
    pure_func < "$item" > "$logname" &
done
wait
