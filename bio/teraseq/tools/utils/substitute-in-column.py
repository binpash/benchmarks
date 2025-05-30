#!/usr/bin/env python3
#
# Substitute a column in tab separated file using 'pattern' in from a tab separated file. First column in the 'pattern' tab separated file is the original string to be replaced, second is the new string.
# Useful to convert, for example, UCSC to Ensembl chromosome names in bed files (primary design). A great resource for chromosome replacement is Devon Ryan's github repository https://github.com/dpryan79/ChromosomeMappings
#

import argparse
import sys

parser = argparse.ArgumentParser(description='Substitute strings in a column with patterns in a substitution table.')
parser.add_argument("-i", "--input", type=str,
                    help="Input file where we'll substitute. Default: stdin")
parser.add_argument("-o", "--output", type=str,
                 help="Output file with finished substitutions. Default: stdout")
parser.add_argument("-t", "--table", type=str,
                 help="Table with two columns for substitution. First column is what to substitute, second what to substitute with. Tab separated.")
parser.add_argument("-s", "--sep", type=str, default="\t",
                help="Column separator in the input file (where we'll substitute). Default: tab")
parser.add_argument("-c", "--column", type=int, default=1,
                help="Column number where strings are to be substituted (where we'll substitute). One-based. Default: 1")

args = parser.parse_args()

## Variables
intab = args.input
otab = args.output
repl_tab = args.table
repl_col = args.column
separ = args.sep

# intab = "in.bed"
# otab = "out.bed"
# repl_tab = "test.tab"
# repl_col = 1
# separ = "\t"

## Input & Output
# Read input
if intab:
    f = open(intab, "r")
else:
    f = sys.stdin

# Get output
if otab:
    fout = open(otab, "w")
else:
    fout = sys.stdout

# Read substitution table and make a dict
subs = open(repl_tab, "r")
subs_lines = subs.readlines()
subs_lines_dict = {}
for line in subs_lines:
    line = line.strip().split("\t")
    if len(line) == 2:
        subs_lines_dict[line[0]] = line[1]

# Read input file lines and replace specified column if found in substitution file
for line in f:
    line = line.strip().split(separ)
    repl = line[repl_col-1]
    if repl in subs_lines_dict.keys(): # If column in substitution, replace and write
        line[repl_col-1] = subs_lines_dict[repl]
        line = separ.join(line)
        fout.write(line + "\n")
    else: # If not found, write as it is
        line = separ.join(line)
        fout.write(line + "\n")

subs.close()

if f is not sys.stdin:
    f.close()

if fout is not sys.stdout:
    fout.close()
