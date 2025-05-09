#!/bin/bash

# These need to be set up:
IN="/srv2/bio4"
IN_NAME="large"
OUT="/srv2/bio4/output"

echo "$IN $IN_NAME $OUT"
echo foo
mkdir -p "$OUT"

# Processing input file
cat "${IN}/${IN_NAME}" | while read s_line; do
    pop="${s_line%%-*}"
    temp="${s_line#*-}"
    sample="${temp%%-*}"

    echo "Processing Sample ${IN}/input/$sample"
    # Correcting labeling of chromosomes so that all are 1,2,3.. instead of chr1,chr2 or chromosome1 etc
    # uniform the chromosomes in the file due to inconsistencies
    samtools view -H "${IN}/input/$sample.bam" | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' \
        | samtools reheader - "${IN}/input/$sample.bam" > "${OUT}/$sample"_corrected.bam

    # Create bai file 
    samtools index -b "${OUT}/$sample"_corrected.bam

    ### Isolating each relevant chromosome based on Gene_locs
    cut -f 2 ./Gene_locs.txt | sort | uniq | while read chr; do
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
