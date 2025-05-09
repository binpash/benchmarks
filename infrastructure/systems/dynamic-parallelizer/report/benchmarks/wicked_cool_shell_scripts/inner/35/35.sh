#!/bin/bash
# fquota--Disk quota analysis tool for Unix; assumes all user
# accounts are >= UID 100
#MAXDISKUSAGE=20000 # In megabytes
MAXDISKUSAGE=2 # In megabytes
for name in $(cut -d: -f1,3 /etc/passwd | awk -F: '$2 > 99 {print $1}')
do
    /bin/echo -n "User $name exceeds disk quota. Disk usage is: "
    # You might need to modify the following list of directories to match
    # the layout of your disk. The most likely change is from /Users to /home.
    find /Users /users /home -xdev -user $name -type f -ls | \
    awk '{ sum += $7 } END { print sum / (1024*1024) " Mbytes" }'
done | awk "\$9 > $MAXDISKUSAGE { print \$0 }"
