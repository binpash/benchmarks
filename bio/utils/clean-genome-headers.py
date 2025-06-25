#!/usr/bin/env python3

import argparse
import sys
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

parser = argparse.ArgumentParser(
    description='Remove content after whitespace in FASTA headers')
parser.add_argument(
    '--fasta', required=True, help='FASTA file')
args = parser.parse_args()

if args.fasta == '-':
    f = sys.stdin
else:
    f = open(args.fasta, 'r')

for l in f:
    if l[0] == ">":
        fields = l.split()
        print(fields[0])
        continue
    print(l, end='')
f.close()

