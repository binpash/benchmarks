#!/usr/bin/env bash
#
# vim:ts=5:sw=5:expandtab
# we have a spaces softtab, that ensures readability with other editors too

# testssl.sh is a program for spotting weak SSL/TLS encryption, ciphers, protocols and some
# vulnerabilities or features. It may or may be not distributed by your distribution.
# The upstream versions are available (please leave the links intact):
#
# Development version       https://github.com/drwetter/testssl.sh
# Stable version            https://testssl.sh
# File bugs at GitHub       https://github.com/drwetter/testssl.sh/issues
#
# Project lead and initiator: Dirk Wetter, copyleft: 2007-today, contributions so far see CREDITS.md
# Main contributions from David Cooper
# Project lead and initiator: Dirk Wetter, copyleft: 2007-today.
# Main contributions from David Cooper. Further contributors see CREDITS.md .
#
# License: GPLv2, see https://opensource.org/licenses/gpl-2.0.php and
# accompanying license "LICENSE.txt". Redistribution + modification under this
# license permitted.
# If you enclose this program or parts of it in your software, it has to be
# accompanied by the same license (see link). Do not violate the license.
# If you do not agree to these terms, do not use it in the first place!
#
# OpenSSL, which is being used and maybe distributed via one of this projects'
# web sites, is subject to their licensing: https://www.openssl.org/source/license.txt
#
# The client simulation data comes from SSLlabs and is licensed to the 'Qualys SSL Labs
# Terms of Use' (v2.2), see https://www.ssllabs.com/downloads/Qualys_SSL_Labs_Terms_of_Use.pdf,
# stating a CC BY 3.0 US license: https://creativecommons.org/licenses/by/3.0/us/
#
# Please note:  USAGE WITHOUT ANY WARRANTY, THE SOFTWARE IS PROVIDED "AS IS".
# USE IT AT your OWN RISK!
# Seriously! The threat is you run this code on your computer and untrusted input e.g.
# could be supplied from a server you are querying.
#
# HISTORY:
# Back in 2006 it all started with a few openssl commands...
# That's because openssl is a such a good swiss army knife (see e.g.
# https://wiki.openssl.org/index.php/Command_Line_Utilities) that it was difficult to resist
# wrapping some shell commands around it, which I used for my pen tests. This is how
# everything started.
# Now it has grown up, it has bash socket support for most features, which has been basically
# replacing more and more functions of OpenSSL and some sockets functions serve as some kind
# of central functions.
#
# WHY BASH?
# Cross-platform is one of the three main goals of this script. Second: Ease of installation.
# No compiling, install gems, go to CPAN, use pip etc. Third: Easy to use and to interpret
# the results.
# /bin/bash including the builtin sockets fulfill all that.  The socket checks in bash may sound
# cool and unique -- they are -- but probably you can achieve e.g. the same result with my favorite
# interactive shell: zsh (zmodload zsh/net/socket -- checkout zsh/net/tcp) too! Oh, and btw.
# ksh93 has socket support too.
# Also bash is quite powerful if you use it appropriately: It can operate on patterns, process lines
# and deal perfectly with regular expressions -- without external binaries.
# /bin/bash though is way more often used within Linux and it's perfect for cross platform support.
# MacOS X has it and also under Windows the MSYS2 extension or Cygwin as well as Bash on Windows (WSL)
# has /bin/bash.
#
# Q: So what's the difference to www.ssllabs.com/ssltest/ or sslcheck.globalsign.com/ ?
# A: As of now ssllabs only check 1) webservers 2) on standard ports, 3) reachable from the
#    internet. And those examples above 4) are 3rd parties. If these restrictions are all fine
#    with you and you need a management compatible rating -- go ahead and use those.
#
# But also if your fine with those restrictions: testssl.sh is meant as a tool in your hand
# and it's way more flexible.  Oh, and did I mention testssl.sh is open source?
#
#################### Stop talking, action now ####################


########### Definition of error codes
#
declare -r ERR_BASH=255            # Bash version incorrect
declare -r ERR_CMDLINE=254         # Cmd line couldn't be parsed
declare -r ERR_FCREATE=253         # Output file couldn't be created
declare -r ERR_FNAMEPARSE=252      # Input file couldn't be parsed
declare -r ERR_NOSUPPORT=251       # Feature requested is not supported
declare -r ERR_OSSLBIN=250         # Problem with OpenSSL binary
declare -r ERR_DNSBIN=249          # Problem with DNS lookup binaries
declare -r ERR_OTHERCLIENT=248     # Other client problem
declare -r ERR_DNSLOOKUP=247       # Problem with resolving IP addresses or names
declare -r ERR_CONNECT=246         # Connectivity problem
declare -r ERR_CLUELESS=245        # Weird state, either though user options or testssl.sh
declare -r ERR_RESOURCE=244        # Resources testssl.sh needs couldn't be read
declare -r ERR_CHILD=242           # Child received a signal from master
declare -r ALLOK=0                 # All is fine


[ -z "${BASH_VERSINFO[0]}" ] && printf "\n\033[1;35m Please make sure you're using \"bash\"! Bye...\033[m\n\n" >&2 && exit $ERR_BASH
[ $(kill -l | grep -c SIG) -eq 0 ] && printf "\n\033[1;35m Please make sure you're calling me without leading \"sh\"! Bye...\033[m\n\n"  >&2 && exit $ERR_BASH
[ ${BASH_VERSINFO[0]} -lt 3 ] && printf "\n\033[1;35m Minimum requirement is bash 3.2. You have $BASH_VERSION \033[m\n\n"  >&2 && exit $ERR_BASH
[ ${BASH_VERSINFO[0]} -le 3 ] && [ ${BASH_VERSINFO[1]} -le 1 ] && printf "\n\033[1;35m Minimum requirement is bash 3.2. You have $BASH_VERSION \033[m\n\n"  >&2 && exit $ERR_BASH

########### Debugging helpers + profiling
#
declare -r PS4='|${LINENO}> \011${FUNCNAME[0]:+${FUNCNAME[0]}(): }'
DEBUGTIME=${DEBUGTIME:-false}                     # https://stackoverflow.com/questions/5014823/how-to-profile-a-bash-shell-script-slow-startup#20855353
DEBUG_ALLINONE=${DEBUG_ALLINONE:-false}           # true: do debugging in one screen (old behavior for testssl.sh and bash3's default
                                                  # false: needed for performance analysis or useful for just having an extra file
DEBUG_ALLINONE=${SETX:-false}                     # SETX as a shortcut for old style debugging, overriding DEBUG_ALLINONE
if [[ "$SHELLOPTS" =~ xtrace ]]; then
     if "$DEBUGTIME"; then
          # separate debugging, doesn't mess up the screen, $DEBUGTIME determines whether we also do performance analysis
          exec 42>&2 2> >(tee /tmp/testssl-$$.log | sed -u 's/^.*$/now/' | date -f - +%s.%N >/tmp/testssl-$$.time)
          # BASH_XTRACEFD=42
     else
          if ! "$DEBUG_ALLINONE"; then
               exec 42>| /tmp/testssl-$$.log
               BASH_XTRACEFD=42
          fi
     fi
fi