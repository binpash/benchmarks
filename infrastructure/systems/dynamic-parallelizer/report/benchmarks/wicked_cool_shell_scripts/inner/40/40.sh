#!/bin/bash
# addusers--Adds new users to the system, including building their
# home directories, copying in default config data, etc.
# For a standard Unix/Linux system, not OS X.
pwfile="/etc/passwd"
shadowfile="/etc/shadow"
gfile="/etc/group"
hdir="/home"

if [ "$(id -un)" != "root" ] ; then
    echo "Error: You must be root to run this command." >&2
    exit 1
fi

echo "Add new user accounts to $(hostname)"
# read -p "Enter the file path containing usernames: " userfile
userfile="$1"

if [ ! -f "$userfile" ]; then
    echo "Error: File '$userfile' not found." >&2
    exit 1
fi

for login in $(cat "$userfile"); do
    uid="$(awk -F: '{ if (big < $3 && $3 < 5000) big=$3 } END { print big + 1 }' "$pwfile")"
    homedir="$hdir/$login"
    gid="$uid"

    # read -p "Enter the full name for $login: " fullname
    # read -p "Enter the shell for $login: " shell
    fullname="New User"
    shell="/bin/bash"

    echo "Setting up account $login for $fullname..."
    echo "${login}:x:${uid}:${gid}:${fullname}:${homedir}:$shell" >> "$pwfile"
    echo "${login}:*:11647:0:99999:7:::" >> "$shadowfile"
    echo "${login}:x:${gid}:$login" >> "$gfile"

    mkdir "$homedir"
    cp -R /etc/skel/.[a-zA-Z]* "$homedir"
    chmod 755 "$homedir"
    chown -R "${login}:${login}" "$homedir"

    # Setting an initial password
    echo "$login:${login}123456789" > tmp_password
    sudo chpasswd < tmp_password
    rm tmp_password
done
