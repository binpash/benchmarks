#!/bin/bash
# encrypt all files in a directory 

mkdir -p "$2"

pure_func() {
    openssl enc -aes-256-cbc -pbkdf2 -iter 20000 -k 'key' -S 1234567890abcdef
}
export -f pure_func

for item in "$1"/*.pcapng; do
    output_name="$2/$(basename "$item").enc"
    pure_func < "$item" > "$output_name" &
done
wait
