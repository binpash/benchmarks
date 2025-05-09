#!/bin/sh
# tag: trigram_rec

pure_func() {
    input=$1
    TEMPDIR=$(mktemp -d)
    tr -sc '[A-Z][a-z]' '[\012*]' > ${TEMPDIR}/${input}.words
    tail +2 ${TEMPDIR}/${input}.words > ${TEMPDIR}/${input}.nextwords
    tail +3 ${TEMPDIR}/${input}.words > ${TEMPDIR}/${input}.nextwords2
    paste ${TEMPDIR}/${input}.words ${TEMPDIR}/${input}.nextwords ${TEMPDIR}/${input}.nextwords2 | sort | uniq -c
    rm -rf ${TEMPDIR}
}

generate_unique_file() {
    local dir="$OUTPUT_DIR"
    local prefix="strace_log"
    local counter_file="$dir/${prefix}"
    if [ ! -f "$counter_file" ]; then
        echo 0 > "$counter_file"
    fi
    local counter
    counter=$(<"$counter_file")
    counter=$((counter + 1))
    echo "$counter" > "$counter_file"
    local filename="$dir/${prefix}_${counter}"
    echo "$filename"
}

export STRACE="strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i"


for input in $(ls $INPUT_FILE | xargs -I arg1 basename arg1)
do
    logfile=$(generate_unique_file)
    strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i cat $INPUT_FILE/$input | grep 'the land of' | pure_func $input | sort -nr | sed 5q > $OUTPUT_DIR/$input.0.out
    
    logfile=$(generate_unique_file)
    strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i cat $INPUT_FILE/$input | grep 'And he said' | pure_func $input | sort -nr | sed 5q > $OUTPUT_DIR/$input.1.out
done
