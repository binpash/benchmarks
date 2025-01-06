#!/bin/bash
# encrypt all files in a directory 
# mkdir -p $2

# pure_func() {
#     openssl enc -aes-256-cbc -pbkdf2 -iter 20000 -k 'key' -S 1234567890abcdef
# }
# export -f pure_func

# for item in $1/*.pcapng;
# do
#     output_name="$2/$(basename $item).enc"
#     cat $item | pure_func > $output_name
# done

# Using GNU parallel:

mkdir -p "$2"

encrypt_file() {
  input_file="$1"
  output_file="$2/$(basename "$input_file").enc"
  openssl enc -aes-256-cbc -pbkdf2 -iter 20000 -k "key" -S "1234567890abcdef" -in "$input_file" -out "$output_file"
}

export -f encrypt_file

find "$1" -type f -name "*.pcapng" | parallel --jobs "$(nproc)" encrypt_file {} "$2"
