#!/bin/sh
# verses with 2 or more, 3 or more, exactly 2 instances of light.

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
    strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i cat $INPUT_FILE/$input | grep -c 'light.\*light'                                 > $OUTPUT_DIR/$input.out0
    
    logfile=$(generate_unique_file)
    strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i cat $INPUT_FILE/$input | grep -c 'light.\*light.\*light'                         > $OUTPUT_DIR/$input.out1
    
    logfile=$(generate_unique_file)
    strace -y -f --seccomp-bpf --trace=fork,clone,%file -o $logfile env -i cat $INPUT_FILE/$input | grep 'light.\*light' | grep -vc 'light.\*light.\*light' > $OUTPUT_DIR/$input.out2
done
