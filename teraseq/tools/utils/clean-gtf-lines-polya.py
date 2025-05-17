#!/usr/bin/env python3

import re
import argparse
import sys
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

def print_active_transcript(tr):
    score = 0;
    if 'start_codon' in tr and tr['start_codon'] == True:
        score += 1
    if 'stop_codon' in tr and tr['stop_codon'] == True:
        score += 1
    if score != 0 and score != 2:
        return

    for l in tr['lines']:
        print(l)

parser = argparse.ArgumentParser(
    description='Keep only specific lines from the GTF. Keep only protein coding transcripts. Only keep lines for "exon", "start_codon" and "stop_codon". Discard transcripts that only have one of start/stop codons defined.')
parser.add_argument(
    '--gtf', required=True, help='Ensembl GTF file')
args = parser.parse_args()

if args.gtf == '-':
    f = sys.stdin
else:
    f = open(args.gtf, 'r')

active_transcript = {'id': ""}
for l in f:
    l = l.rstrip()
    if l[0] == "#":
        print(l)
        continue

    fields = l.split("\t")
    if fields[2] not in ["exon", "start_codon", "stop_codon"]:
        continue

    # Make sure that the lines are of the following type.
    discard = True
    ps = []
    ps.append(re.compile('transcript_biotype "protein_coding"')) # original by manolis
    # Append all long non-coding RNAs which might have polyA as well and were found in cip.decap sample
    ps.append(re.compile('transcript_biotype "lincRNA"')) # original by manolis
#    ps.append(re.compile('transcript_biotype "non_stop_decay"'))
    ps.append(re.compile('transcript_biotype "antisense_RNA"'))
    ps.append(re.compile('transcript_biotype "bidirectional_promoter_lncRNA"'))
#    ps.append(re.compile('transcript_biotype "nonsense_mediated_decay"'))
    ps.append(re.compile('transcript_biotype "processed_transcript"'))
    ps.append(re.compile('transcript_biotype "retained_intron"'))
    ps.append(re.compile('transcript_biotype "sense_intronic"'))
    ps.append(re.compile('transcript_biotype "sense_overlapping"'))

    for p in ps:
        if p.search(l):
            discard = False
    if discard:
        continue

    # Make sure that the lines don't contain the following tags.
    ps = []
    ps.append(re.compile('tag "mRNA_start_NF";')) # the mRNA start could not be confirmed.
    ps.append(re.compile('tag "mRNA_end_NF";')) # the mRNA end could not be confirmed.
    ps.append(re.compile('tag "cds_start_NF";')) # the coding region start could not be confirmed.
    ps.append(re.compile('tag "cds_end_NF";')) # the coding region end could not be confirmed.
    for p in ps:
        if p.search(l):
            discard = True
    if discard:
        continue

    # Make sure that all coding transcripts have both start and stop codon.
    tr_p = re.compile('transcript_id "(.+?)"')
    m = tr_p.search(l)
    if m:
        tid = m.group(1)

    if tid != active_transcript['id']:
        if active_transcript['id'] != "":
            print_active_transcript(active_transcript)
        active_transcript = {'id': tid, 'lines': []}

    if fields[2] == "start_codon":
        active_transcript["start_codon"] = True

    if fields[2] == "stop_codon":
        active_transcript["stop_codon"] = True

    active_transcript["lines"].append(l)

f.close()

if active_transcript['id'] != "":
    print_active_transcript(active_transcript)
