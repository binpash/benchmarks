#!/usr/bin/env python3
#
# # Add new SAM flag with number of mappings per read name
# Add number of mappings per read by name. The goal is to avoid dependance on mapper MAPQ (or other way it assigns unique mappings). A very beta version. So far cannot handle chimeric and paired-end reads. 
#
#TODO: Add support for paired-end reads (paired, singletons). 
#TODO: Add support for chimeric reads. 
#

import argparse
import sys
from collections import defaultdict
from collections import Counter
import re

# Function to see head of an object
#from itertools import islice
#
#def take(n, iterable):
#    "Return first n items of the iterable as a list"
#    return list(islice(iterable, n))

# Function to test SAM flag
def checksamflag(checkflag, samflag):
    "Return True of False is checked flag is in sam flag"
    return str(format(int(checkflag), "b"))[::-1].find(str(1)) in [m.start() for m in re.finditer(str(1), str(format(int(samflag), "b")[::-1]))]

parser = argparse.ArgumentParser(description='Count number of mappings per read.')
parser.add_argument("-i", "--input", 
					help="Input SAM file. Default: stdin")
parser.add_argument("-o", "--output", 
				 help="Output SAM file. Default: stdout")
parser.add_argument("-t", "--tag", type=str, default="X0", 
				help="New tag with number of mappings to be added. Default: X0.")
parser.add_argument("-w", "--overwrite", action='store_true', default = False,  
				help="Do you want to overwrite if --tag already exists? Default: Don't overwrite.")

args = parser.parse_args()

## Variables
input = args.input
#input = "test.sam"

output = args.output
#output = "out.sam"

tag_add = args.tag # "X0"
#tag_add = "X0"

# Do you want to overwrite if tag already exists? Default: False
over = args.overwrite 
#over = "True"

## Input & Output
# Read input
if input:
    f=open(input, "r").read()
else:
    f=sys.stdin.read()

lines=f.split('\n')

# Get output
if output:
    fout=open(output, "w")
else:
    fout=sys.stdout

index = defaultdict(list)
counts = {}   

for line in range (len(lines)):
#    print(lines[line])
    if lines[line].startswith('@'): # Header save
#        print("It's a header!")
        fout.write(lines[line] + '\n')
    else:
        if not lines[line].strip(): # Skip empty lines
            continue
        elif checksamflag(4, lines[line].split('\t')[1]) or checksamflag(2048, lines[line].split('\t')[1]) or checksamflag(1, lines[line].split('\t')[1]): # Skip if read is unmapped, chimeric or paired
#            print("Read unmapped, chimeric or paired, skipping.")
            fout.write(lines[line] + '\n')
        else:
            if [i for i in lines[line].split('\t')[11:] if i.startswith(tag_add+":")]:
                if over:
#                    print("Warning! Found already existing " + tag_add + " and \'overwrite\' is set to: " + over + ". We are going to overwrite.")
                    index[lines[line].split('\t')[0]].append(line) # Reads
                else:
                    sys.exit("Found already existing " + tag_add + " but 'overwrite' is set to False, not going to overwrite and exiting. Set '--overwrite' to enable overwriting or choose a different flag to be added.")
            else:
#                print("Haven't found " + tag_add + " in the SAM line, going to add it.")
                index[lines[line].split('\t')[0]].append(line) # Reads

# Simple counting of number of values for each key
for key, value in index.items():
    counts[key] = len(value)

# Test to see what's at least twice
#for key, value in counts.items():
#    if value == 2:
#        print(key)
#        print(counts[key])

for key, value in counts.items():   

    lines_read = [lines[i] for i in index[key]]
    
    for element in lines_read:
#        print(element)
        z = re.search('\t'+tag_add+':i:[0-9]+', element)
        if z is None: # If there is not match
            fout.write(element + '\t' + tag_add + ':i:' + str(value) + '\n') # Add "number of alignments" tag
        else: # If there is a match
            # Get everything before match & after match and replace the value in between 
            fout.write(element[:z.span()[0] + len(tag_add) + 4] + str(value) + element[z.span()[1] + 1:] + '\n')        

if output:
    fout.close()
