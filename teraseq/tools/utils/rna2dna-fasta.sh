#!/bin/bash
#
# converts RNA fasta to DNA fasta 
# converts u/U to t/T on every which doesn't start with ">"
#

if [[ $# -eq 0 ]] ; then
    echo 'Please, enter exactly one argument with RNA (Uu) fasta to convert.'
    exit 0
fi

if [[ $# -eq 1 ]] ; then
#	echo 'Converting ${1} to DNA (U->T or u->t)'
	sed '/^[^>]/ y/uU/tT/' $1 
else
	echo 'Please enter EXACTLY one argument - one RNA (Uu) fasta to be converted to DNA (Tt).'
fi
