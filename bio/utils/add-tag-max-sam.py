#!/usr/bin/env python3
#
# Get the max tag value for a single read and add to all the corresponding reads
# Initially designed for adding count of primary alignments to XP:i tag for Minimap2 alignments based on read name and a value of ms:i tag
#
# Input is sam file, output is sam file (unsorted)
#
# It's quite fast but consumes a lot of RAM. For 2.5 GB BAM it needs ~13 GB RAM 
#

import argparse
import sys
import time
from collections import defaultdict

parser = argparse.ArgumentParser(description='Get a count of the maximum value of a specified SAM tag and output the count to a new tag.')
parser.add_argument("-i", "--input", 
					help="Input SAM file. Default: stdin")
parser.add_argument("-o", "--output", 
				 help="Output SAM file. Default: stdout")
parser.add_argument("-t", "--tag", type=str, default="ms", 
				help="Tag name to get the maximum from. Default: ms.")
parser.add_argument("-n", "--newtag", type=str, default="XP", 
                help="Tag name to determine 'maximum value hits' (~best hits) from --tag. Can be 0 for no value, 1 for the max value hits, 2 for all the other hits. Default: XP.")
parser.add_argument("-n2", "--newtag2", type=str, default="XN", 
				help="Tag name to add the count of maximum value of --tag. Can be 0 for no value, any number for number of max values. Default: XN.")

args = parser.parse_args()

tag_scan = args.tag # "ms"
tag_add = args.newtag # "XP"
tag_add2 = args.newtag2 # "XN"

t0 = time.time()

# Read input
if args.input:
    f=open(args.input, "r").read()
else:
    f=sys.stdin.read()

lines=f.split('\n')

# Get output
if args.output:
    fout=open(args.output, "w")
else:
    fout=sys.stdout

reads = defaultdict(list)
values = defaultdict(list)
read_notag = []
read_singletag = []
read_multitag = []

# Get line number (index) of reads w/wo tag and for those with a tag get the values
for line in range (len(lines)):
    if lines[line].startswith('@'):
 #       print("It's a header!")
        fout.write(lines[line] + '\n')
    elif [i for i in lines[line].split('\t')[11:] if i.startswith(tag_scan+":")]:
#        print("Found it!")
        reads[lines[line].split('\t')[0]].append(line) # Get index of reads
        values[lines[line].split('\t')[0]].append([i for i in lines[line].split('\t')[11:] if i.startswith(tag_scan+":")][0].rsplit(':', 1)[1]) # Get read tag values
    else:
#        print("Didn't find it")
        if not lines[line].strip(): # Skip empty lines
            continue
        else:
            read_notag.append(line) # Get index of reads with no tag

reads_l = [v for k,v in reads.items()] # get only values from the dictinary
values_l = [v for k,v in values.items()] # get only values from the dictinary

read_multitag = [i for i, x in enumerate([len(x) for x in reads_l]) if x > 1] # Get read with tags more than once https://stackoverflow.com/questions/6294179/how-to-find-all-occurrences-of-an-element-in-a-list
read_singletag = [i for i, x in enumerate([len(x) for x in reads_l]) if x == 1] # Get read with only one tag

for x in read_multitag:
    values_x = [int(i) for i in values_l[x]] # Make sure the values are int and not str otherwise max() doesn't work properly
    max_val = values_x.count(max(values_x))
    max_val_ind = [i for i, x in enumerate(values_x) if x == max(values_x)] # Get possition of all the max values

    # Get rows by indeces and add first tag, for the rest don't - index to row number (different!) - index is from sublist, row number from total read lines
    reads_l_max = [reads_l[x][i] for i in max_val_ind] # Get reads with the "best hit"
    reads_l_rest = [reads_l[x][i] for i in list(range(0, len(values_x))) if i not in max_val_ind] # Get the other read

    for i in reads_l_max:
        fout.write(lines[i] + '\t' + tag_add + ':i:1' + '\t' + tag_add2 + ':i:' + str(max_val) + '\n') # Add "best hit" tag
    for i in reads_l_rest:
        fout.write(lines[i] + '\t' + tag_add + ':i:2' + '\t' + tag_add2 + ':i:' + str(max_val) + '\n') # Add "best hit" tag

flat_list_single = [item for sublist in [reads_l[i] for i in read_singletag] for item in sublist] # Unlist the stupid list https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists

for i in flat_list_single:
	fout.write(lines[i] + '\t' + tag_add + ':i:1' + '\t' + tag_add2 + ':i:1' + '\n') # Add "best hit" tag
for i in read_notag:
	fout.write(lines[i] + '\t' + tag_add + ':i:0' + '\t' + tag_add2 + ':i:0' + '\n') # Add "no value" tag

if args.output:
    fout.close()

#print("Don't forget to resort the SAM!")
