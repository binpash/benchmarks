#!/usr/bin/env python3

"""
Reads a FASTQ file and checks if the 5' end of the sequenced can be located
with cutadapt.
"""

# TODO: Group reads by reference, trim the same adapter from the group of reads, move to the next one

import argparse
import sys
import subprocess
import re
import os
import pysam
import numpy as np
from numpy import random
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
matplotlib.use('pdf')


def get_file_object(filename):
    """
    Returns a file objects that reads from filename; reads from STDIN if
    filename is -
    """
    if filename == "-":
        return sys.stdin
    return open(filename, "r")


def write_fasta_entry_to_file(filename, name, seq):
    """
    Creates a file and writes the fasta entry in there.
    """
    w = open(filename, "w")
    w.write(">" + name + "\n")
    w.write(seq)
    w.close()
    return filename


def reads_with_adapter_from_cutadapt(filename, adaptor):
    """
    Run cutadapt on the filename and returns the number of reads where the
    adaptor was found in
    """
    # Prepare the cutadapt options.
    cutadapt_args = ['-g', adaptor,
                     '--overlap', '31',
                     '--minimum-length', '25',
                     '--error-rate', '0.29',
                     '--output', 'foo.txt',
                     '--discard-untrimmed',
                     filename]

    # Run cutadapt to see if we recognise the adaptor in the sequence.
    res = subprocess.run(['cutadapt'] + cutadapt_args,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL)

    m = re.search(r'Reads\swith\sadapters:\s+(\d+)', str(res.stdout))
    if m is None:
        print(res.stdout, file=sys.stderr)
        exit(1)

    os.remove("foo.txt") # Clean

    return int(m.group(1))


# Read command line parameters
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-b", "--bam", required=True,
                    help="Input BAM file")
parser.add_argument("-f", "--fasta", required=True,
                    help="Input FASTA with reference sequences")
parser.add_argument("-q", "--fastq", required=True,
                    help="Input FASTQ with previously adaptor trimmed reads")
parser.add_argument("-o", "--figfile", required=True,
                    help="Output figure file")
args = parser.parse_args()

# Define global vars.
ref_seqs = {}     # Holds names and sequences of reference sequences.
reads_with_adapter = {}  # Holds the names of previously trimmed reads.
found = {}    # Number of reads with identified adaptor.
scanned = {}  # Number of reads scanned.
reads_with_adapter_found = {}   # Number of reads with identified adaptor from previously trimmed reads.
reads_with_adapter_scanned = {}  # Number of reads scanned from previously trimmed reads.

# Read FASTA file with reference sequences.
fa = get_file_object(args.fasta)
for ln in fa:
    if ln[0] != '>':
        continue
    name = ln[1:].strip()
    name = name.split()[0] # Split name after first space
    seq = next(fa).strip()
    ref_seqs[name] = seq.upper()

# Read FASTQ file with reads where adaptor was previously identified.
fq = get_file_object(args.fastq)
for ln in fq:
    if ln[0] != '@':
        continue
    name = ln[1:].strip()
    reads_with_adapter[name] = True
    next(fq)
    next(fq)
    next(fq)

if not os.path.exists("tmp"):
    os.mkdir("tmp")

# Loop on the BAM file
samfile = pysam.AlignmentFile(args.bam, "rb")
for read in samfile:
    start = read.reference_start
    query_name = read.query_name
    query_seq = read.query_sequence
    ref_seq = ref_seqs.get(read.reference_name)
    if ref_seq is None:
        print("Error: No reference seq for: " + read.reference_name,
              file=sys.stderr)
        exit(1)

    # Skip very short reads (100 bp)
    min_len = 100
    if len(query_seq) < min_len:
        continue

    # Only keep reads that align close to the TSS (first 40 bp)
    dist = 40
    if start > dist:
        continue

    # Assume adaptor is a sequence at the 5' end of the reference sequence.
    adaptor_len = 58
    adaptor = ref_seq[0:adaptor_len]
    # Note: We could anchor the adapter to the 5' end but since there is a possibility of having some misalignments at the 5' we won't do it now
    #adaptor = "X" + adaptor

    # Write the read sequence to a temporary file.
    temp_file = "tmp/" + query_name + ".fasta"
    write_fasta_entry_to_file(temp_file, query_name, query_seq)

    # Run cutadapt and get the number of reads with the adaptor.
    # NOTE: Comment the line below and uncomment the next for quickly testing
    # the script.
    cnt_reads_with_ad = reads_with_adapter_from_cutadapt(temp_file, adaptor)
    # cnt_reads_with_ad = random.randint(0, 2)
    os.remove(temp_file)

#    base = 5
    base = 1
    start = base * round(start/base)  # Round to the nearest multiple of 'base'.
    # Initialize dictionaries.
    if start not in scanned:
        scanned[start] = 0
        found[start] = 0
        reads_with_adapter_scanned[start] = 0
        reads_with_adapter_found[start] = 0

    if reads_with_adapter.get(query_name, False):
        reads_with_adapter_scanned[start] += 1
        reads_with_adapter_found[start] += cnt_reads_with_ad
    else:
        scanned[start] += 1
        found[start] += cnt_reads_with_ad

#    if sum(scanned.values()) % 10000 == 0:
#        print("Processed: " + str(sum(scanned.values())), file=sys.stderr)

    # if sum(scanned.values()) > 0 and sum(scanned.values()) % 200000 == 0:
    #     break
samfile.close()

minkey = min(scanned.keys())
maxkey = max(scanned.keys())
print("\t".join(['start', 'scanned', 'found',
                 'reads_with_adapter_scanned', 'reads_with_adapter_found']))
for k in range(minkey, maxkey+1):
    print("\t".join([str(k),
                     str(scanned.get(k, 0)),
                     str(found.get(k, 0)),
                     str(reads_with_adapter_scanned.get(k, 0)),
                     str(reads_with_adapter_found.get(k, 0))]))

s = np.array([scanned.get(k, 0) for k in range(minkey, maxkey+1)])
f = np.array([found.get(k, 0) for k in range(minkey, maxkey+1)])
ss = np.array([reads_with_adapter_scanned.get(k, 0) for k in range(minkey, maxkey+1)])
ff = np.array([reads_with_adapter_found.get(k, 0) for k in range(minkey, maxkey+1)])
df = pd.DataFrame(data={
    'wo_rel5': f/s,
    'w_rel5': ff/ss}).dropna()
sns_plt = sns.lineplot(data=df)
sns_plt.set(ylim=(0, 1))
sns_plt.set(xlabel='Alignment_start', ylabel='Sensitivity')
plt.savefig(args.figfile)
