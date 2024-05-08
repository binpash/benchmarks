#!/bin/bash
# First command is almost always a generator

grep -iv ': starting\|kernel: .*: Power Button\|watching system buttons\|Stopped Cleaning Up\|Started Crash recovery kernel' /var/log/messages /var/log/syslog /var/log/* 2> /dev/null | 
    grep --regex 'recover[a-z]*\|power[a-z]*\|shut[a-z ]*down\|rsyslogd\|ups' > ${OUT}/shutdown.log 
