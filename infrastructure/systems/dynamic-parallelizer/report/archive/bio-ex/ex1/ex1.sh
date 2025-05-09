#!/bin/bash

# Viewing the first few lines of chr20.bed
echo "First lines of chr20.bed:"
head chr20.bed

# Count lines in chr20.bed and chr20_flank500.bed
echo "Line counts:"
wc -l chr20.bed chr20_flank500.bed

# Sort ATF3_chr20.bed file
echo "Sorting ATF3_chr20.bed..."
cat ATF3_chr20.bed | bedtools sort -i > ATF3_chr20_sorted.bed

# Intersect ATF3 binding sites with promoter regions for chromosome 20
echo "Finding common regions for Chromosome 20..."
bedtools intersect -wa -a ATF3_chr20_sorted.bed -b chr20_flank500.bed > common_chr20.bed

# Repeat sorting and intersecting for chromosome 19
echo "Sorting ATF3_chr19.bed..."
bedtools sort -i ATF3_chr19.bed > ATF3_chr19_sorted.bed

echo "Finding common regions for Chromosome 19..."
bedtools intersect -wa -a ATF3_chr19_sorted.bed -b chr19.bed > common_chr19.bed

echo "Analysis completed. Check the Exercise1 directory for results."
