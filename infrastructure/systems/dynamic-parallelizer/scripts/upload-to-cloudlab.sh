#!/bin/bash

ip=${1?IP of remote machine not given}
user=${2?user in remote machine}

## Optionally the caller can give us a private key for the ssh
key=$3
if [ -z "$key" ]; then
    key_flag=""
else
    key_flag="-i ${key}"
fi

export pash_spec_dir=$(git rev-parse --show-toplevel --show-superproject-working-tree)

excludes="--exclude .git/modules/deps --exclude *.tar.gz"
## Upload the whole directory to a cloudlab machine
rsync --rsh="ssh -p 22 ${key_flag}" --progress ${excludes} -p -r "${pash_spec_dir}" "${user}@${ip}:/users/${user}"
    


