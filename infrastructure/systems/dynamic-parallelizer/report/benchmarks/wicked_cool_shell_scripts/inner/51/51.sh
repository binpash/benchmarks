#!/bin/bash
# backup--Creates either a full or incremental backup of a set of defined
# directories on the system. By default, the output file is compressed and
# saved in /tmp with a timestamped filename. Otherwise, specify an output
# device (another disk, a removable storage device, or whatever else floats
# your boat).

compress="bzip2" # Change to your favorite compression app.
# inclist="/tmp/backup.inclist.$(date +%d%m%y)" # Changed this to custom file
inclist="$1/backup.inclist.$$" # Change this to a local temp file.
# output="/tmp/backup.$(date +%d%m%y).bz2" # Changed this to custom file
output="$1/backup.$(date +%d%m%y).bz2"
# tsfile="$HOME/.backup.timestamp" # Changed this to custom file
tsfile="$2"
# btype="incremental" # Default to an incremental backup.
btype="full" # Default to a full backup.

# trap "/bin/rm -f $inclist" EXIT

########## Main code section begins here ###########

if [ "$1" = "-f" ]; then
    btype="full"
fi

echo "Doing $btype backup, saving output to $output"

timestamp="$(date +'%m%d%I%M')" # Grab month, day, hour, minute from date.
# Curious about date formats? "man strftime"

if [ "$btype" = "incremental" ] ; then
    if [ ! -f $tsfile ] ; then
        echo "Error: can't do an incremental backup: no timestamp file" >&2
        exit 1
    fi

    find "$RESOURCE_DIR/h" -depth -type f -newer $tsfile -user $(id -u) | \
    tar --format=posix -cf - -T - | $compress > $output
    failure="$?"
else
    find "$RESOURCE_DIR/h" -depth -type f -newer $tsfile -user $(id -u) | \
    tar --format=posix -cf - -T - | $compress > $output
    failure="$?"
fi

if [ "$failure" = "0" ] ; then
    touch -t $timestamp $tsfile
fi
