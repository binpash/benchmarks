#!/usr/bin/env python3

import sys
import re

## Author: https://gist.github.com/jqsunac/aeca04ee2c5b5cc53ad795b660edd6c3
##
## Description:
##   Calculation of non-overlapping exon length with GFF file for each gene.
## 
##   e.g) Gene G has four transcripts: Ga, Gb, Gc, and Gd. The four transcripts
##        have different numbers of exons and different combinations of exons.
##        Some regions are shared with the four transcripts, and some regions
##        are only used by a single transcript. These information are saved in
##        GTF file. This script is used for calculating the non-overlapping
##       exon (cDNA) length with GFF file.
## 
##    transcript Ga     ========-----==============-------============= 
##    transcript Gb     ========-----==============
##    transcript Gc     ========--------------------------=============
##    transcript Gd     ========--------------==========--=============
##    
##    cDNA              ========     ===================  =============
##
## Usage:
##
##   python calc_cdna_len.py gene_id ath.gtf > cds_length.txt
##
##     gene_id: 'gene_id' can be changed to 'transcript_id', or 'gene_name'
##              according your purpose and your GTF file.
##
##     ath.gtf: File path to the GTF file.
## 


def calc_cdna_len(attr_pattern, gff_file):
    #
    # Calculate the length of non-overlapping exons for each gene.
    #

    idptn = re.compile(attr_pattern + ' "([^"]+)";')
    gene_length = {}

    # save the coordinates of the positions of all exons.
    with open(gff_file, 'r') as gfffh:
        for buf in gfffh:
            if buf[0:1] == '#':
                continue
            buf_records = buf.split('\t')
            if buf_records[2] == 'exon':
                mat = idptn.search(buf_records[8])
                if mat:
                    gene_name = mat.group(1)
                    if gene_name not in gene_length:
                        gene_length[gene_name] = []
                    gene_length[gene_name].append([int(buf_records[3]), int(buf_records[4])])

    # find the maximum and minimum coodinates of the position of exons for each gene.
    gene_length_max = {}
    gene_length_min = {}
    for gene_name in gene_length.keys():
        for exon_range in gene_length[gene_name]:
            # define
            if gene_name not in gene_length_max:
                gene_length_max[gene_name] = max(exon_range)
            if gene_name not in gene_length_min:
                gene_length_min[gene_name] = min(exon_range)
            # update
            if max(exon_range) > gene_length_max[gene_name]:
                gene_length_max[gene_name] = max(exon_range)
            if min(exon_range) < gene_length_min[gene_name]:
                gene_length_min[gene_name] = min(exon_range)
    
    # generate an list contained 0 or 1 for each gene.
    # '0' means that the position did not become a exon region,
    # '1' means that the position became a exon region at least once.
    gene_length_bits = {}
    for gene_name in gene_length.keys():
        if gene_name not in gene_length_bits:
            gene_length_bits[gene_name] = [0] * (gene_length_max[gene_name] - gene_length_min[gene_name] + 1)
        for exon_range in gene_length[gene_name]:
            for i in range(min(exon_range), max(exon_range) + 1):
                pos_i = i - gene_length_min[gene_name]
                gene_length_bits[gene_name][pos_i] = 1
    
    # calculate the total values of the list contained 0 or 1.
    gene_length_final = {}
    for gene_name in gene_length.keys():
        gene_length_final[gene_name] = sum(gene_length_bits[gene_name])

    return gene_length_final


if __name__ == '__main__':
    gene_length = calc_cdna_len(sys.argv[1], sys.argv[2])
    
    for gene_name, gene_len in gene_length.items():
        print(gene_name + '\t' + str(gene_len))
