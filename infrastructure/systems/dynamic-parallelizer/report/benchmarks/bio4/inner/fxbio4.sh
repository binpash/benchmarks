#!/bin/bash

# These need to be set up:
BIODIR="${BIODIR:-$PASH_SPEC_TOP/report/benchmarks/bio4}"
IN="${IN:-$BIODIR/../../resources/bio4}"
IN_NAME=$1
OUT="${OUT:-$BIODIR/output}"

echo "$IN $IN_NAME $OUT"
#cd $IN
mkdir -p "$OUT"

# Processing input file
for line in $(cat "${BIODIR}/${IN_NAME}")
do
	#sample=$(echo "$line" | cut -d " " -f 2)
	#pop=$(echo "$line" | cut -f 1 -d " ")
	#link=$(echo "$line" | cut -f 3 -d " ")
	#pop=$(echo "$line" | sed 's/^\([^-]*\)-.*/\1/')
	#sample=$(echo "$line" | sed 's/^[^-]*-\([^-]*\)-.*/\1/')
	pop="${line%%-*}"
	temp="${line#*-}"
	sample="${temp%%-*}"

	echo "Processing Sample ${IN}/$sample"
	# Correcting labeling of chromosomes so that all are 1,2,3.. instead of chr1,chr2 or chromosome1 etc
	# uniform the chromosomes in the file due to inconsistencies
	samtools view -H "${IN}/$sample.bam" | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' \
		| samtools reheader - "${IN}/$sample.bam" > "${OUT}/$sample"_corrected.bam

    # Create bai file 
    samtools index -b "${OUT}/$sample"_corrected.bam

    ### Isolating each relevant chromosome based on Gene_locs
    for chr in $(cut -f 2 "$BIODIR"/Gene_locs.txt | sort | uniq); do
	    echo "Isolating Chromosome $chr from sample ${OUT}/$sample"
	    samtools view -b "${OUT}/$sample"_corrected.bam chr"$chr" > "${OUT}/$pop"_"$sample"_"$chr".bam
	    echo "Indexing Sample $pop'_'${OUT}/$sample"
	    samtools index -b "${OUT}/$pop"_"$sample"_"$chr".bam
	    # Sleep command can be uncommented if needed
	    # sleep 2
    done

    # Clean up commands can be uncommented if file removal is desired
    # rm "${OUT}/$sample"_corrected.bam
    # rm "${OUT}/$sample"_corrected.bam.bai
    # rm "${OUT}/$sample".bam
done
