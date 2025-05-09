#!/bin/sh
# deleteuser--Deletes user accounts without a trace.
# Not for use with OS X.

homedir="/home"
pwfile="/etc/passwd"
shadow="/etc/shadow"
newpwfile="/etc/passwd.new"
newshadow="/etc/shadow.new"
suspend="$(which suspenduser)"
locker="/etc/passwd.lock"

if [ -z $1 ] ; then
    echo "Usage: $0 file" >&2
    exit 1
elif [ "$(whoami)" != "root" ] ; then
    echo "Error: you must be 'root' to run this command." >&2
    exit 1
fi

file=$1

for username in $(cat "$file"); do
 #   $suspend $username # Suspend their account while we do the dirty work.

    uid="$(grep -E "^${username}:" $pwfile | cut -d: -f3)"

    if [ -z $uid ] ; then
        echo "Error: no account $username found in $pwfile" >&2
        continue
    fi

    # Remove the user from the password and shadow files.
    grep -vE "^${username}:" $pwfile > $newpwfile
    grep -vE "^${username}:" $shadow > $newshadow

    lockcmd="$(which lockfile)" # Find lockfile app in the path.

    if [ ! -z $lockcmd ] ; then # Let's use the system lockfile.
        eval $lockcmd -r 15 $locker
    else # Ulp, let's do it ourselves.
        while [ -e $locker ] ; do
            echo "waiting for the password file"
            sleep 1
        done
        touch $locker # Create a file-based lock.
    fi

    mv $newpwfile $pwfile
    mv $newshadow $shadow
    rm -f $locker # Click! Unlocked again.
    chmod 644 $pwfile
    chmod 400 $shadow

    # Now remove home directory and list anything left.
    rm -rf $homedir/$username

    echo "Files still left to remove (if any):"
    find / -uid $uid -print 2>/dev/null | sed 's/^/ /'
    echo ""
    echo "Account $username (uid $uid) has been deleted, and their home directory "
    echo "($homedir/$username) has been removed."

done

exit 0
