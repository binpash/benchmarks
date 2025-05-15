#!/usr/bin/env python

import re
import argparse
import sys
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

def print_active_transcript(tr):
#    score = 0;
#    if 'start_codon' in tr and tr['start_codon'] == True:
#        score += 1
#    if 'stop_codon' in tr and tr['stop_codon'] == True:
#        score += 1
#    if score != 0 and score != 2:
#        return

    for l in tr['lines']:
        print(l)

parser = argparse.ArgumentParser(
    description='Keep only specific lines from the GTF. Keep only non pseudo genes and non rRNA. Only keep lines for "exon", "start_codon" and "stop_codon". DO NOT discard transcripts that only have one of start/stop codons defined because we have non-coding RNAs which do not have to have the start/stop codong.')
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

    ps.append(re.compile('transcript_biotype "3prime_overlapping_ncRNA"'))
    ps.append(re.compile('transcript_biotype "antisense_RNA"'))
    ps.append(re.compile('transcript_biotype "bidirectional_promoter_lncRNA"'))
    ps.append(re.compile('transcript_biotype "IG_C_gene"'))
    ps.append(re.compile('transcript_biotype "IG_D_gene"'))
    ps.append(re.compile('transcript_biotype "IG_J_gene"'))
    ps.append(re.compile('transcript_biotype "IG_V_gene"'))
    ps.append(re.compile('transcript_biotype "lincRNA"'))
    ps.append(re.compile('transcript_biotype "macro_lncRNA"'))
    ps.append(re.compile('transcript_biotype "miRNA"'))
    ps.append(re.compile('transcript_biotype "misc_RNA"'))
    ps.append(re.compile('transcript_biotype "Mt_tRNA"'))
    ps.append(re.compile('transcript_biotype "non_coding"'))
    ps.append(re.compile('transcript_biotype "nonsense_mediated_decay"'))
    ps.append(re.compile('transcript_biotype "non_stop_decay"'))
    ps.append(re.compile('transcript_biotype "processed_transcript"'))
    ps.append(re.compile('transcript_biotype "protein_coding"'))
    ps.append(re.compile('transcript_biotype "retained_intron"'))
    ps.append(re.compile('transcript_biotype "sense_overlapping"'))
    ps.append(re.compile('transcript_biotype "scaRNA"'))
    ps.append(re.compile('transcript_biotype "scRNA"'))
    ps.append(re.compile('transcript_biotype "sense_intronic"'))
    ps.append(re.compile('transcript_biotype "snoRNA"'))
    ps.append(re.compile('transcript_biotype "snRNA"'))
    ps.append(re.compile('transcript_biotype "sRNA"'))
    ps.append(re.compile('transcript_biotype "TEC"'))
    ps.append(re.compile('transcript_biotype "TR_C_gene"'))
    ps.append(re.compile('transcript_biotype "TR_D_gene"'))
    ps.append(re.compile('transcript_biotype "TR_J_gene"'))
    ps.append(re.compile('transcript_biotype "TR_V_gene"'))
    ps.append(re.compile('transcript_biotype "vaultRNA"'))

    for p in ps:
        if p.search(l):
            discard = False
    if discard:
        continue

    # Make sure that the lines don't contain the following tags.
    # https://www.gencodegenes.org/pages/tags.html
    ps = []
    ps.append(re.compile('tag "mRNA_start_NF";')) # the mRNA start could not be confirmed.
    ps.append(re.compile('tag "mRNA_end_NF";')) # the mRNA end could not be confirmed.
    ps.append(re.compile('tag "cds_start_NF";')) # the coding region start could not be confirmed; this one probably should not be used for total...
    ps.append(re.compile('tag "cds_end_NF";')) # the coding region end could not be confirmed; this one probably should not be used for total...
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
