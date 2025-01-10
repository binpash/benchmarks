#!/bin/bash

# cd "$(dirname "$0")"/.. || exit 1

# IN="$1"
# IN_NAME="$2"
# OUT="$3"

# cat ${IN_NAME}|while read s_line;
#   do
#     sample=$(echo $s_line |cut -d " " -f 2);
#     pop=$(echo $s_line |cut -f 1 -d " ");
#     link=$(echo $s_line |cut -f 3 -d " ");
#     ### correcting labeling of chromosomes so that all are 1,2,3.. instead of chr1,chr2 or chromosome1 etc
#     echo "Processing Sample $sample";
#     # uniform the chromosomes in the file due to inconsistencies
#     samtools view -H "${IN}/$sample".bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' \
#       | samtools reheader - "${IN}/$sample".bam > "${OUT}/$sample"_corrected.bam ;
#     # create bai file 
#     samtools index -b "${OUT}/$sample"_corrected.bam ;
#     ### Isolating each relevant chromosome based on Gen_locs
#     cut -f 2 ./Gene_locs.txt |sort |uniq |while read chr;
#   do  
#     echo 'Isolating Chromosome '$chr' from sample '${OUT}/$sample',  ';
#     samtools view -b "${OUT}/$sample"_corrected.bam chr"$chr" > "${OUT}/$pop"_"$sample"_"$chr".bam ;
#     echo 'Indexing Sample '$pop'_'${OUT}/$sample' ';
#     samtools index -b "${OUT}/$pop"_"$sample"_"$chr".bam;
#   done;
#   #rm "${OUT}/$sample"_corrected.bam;
#   #rm "${OUT}/$sample"_corrected.bam.bai;
#   #rm "${OUT}/$sample".bam
# done;

cd "$(dirname "$0")"/.. || exit 1

IN="$1"
IN_NAME="$2"
OUT="$3"

process_sample() {
    s_line=$1
    sample=$(echo "$s_line" | awk '{print $2}')
    pop=$(echo "$s_line" | awk '{print $1}')
    link=$(echo "$s_line" | awk '{print $3}')

    echo "Processing Sample $sample"

    # Correct labeling of chromosomes and reheader BAM file
    samtools view -H "${IN}/${sample}.bam" | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' \
        | samtools reheader - "${IN}/${sample}.bam" > "${OUT}/${sample}_corrected.bam"

    # Create index file for the corrected BAM
    samtools index -b "${OUT}/${sample}_corrected.bam"

    # Isolate chromosomes based on Gen_locs
    while read -r chr; do
        echo "Isolating Chromosome ${chr} from sample ${OUT}/${sample}"
        samtools view -b "${OUT}/${sample}_corrected.bam" chr"${chr}" > "${OUT}/${pop}_${sample}_${chr}.bam"
    done < <(cut -f 2 ./Gene_locs.txt | sort -u)

    # Index isolated chromosome BAM files
    while read -r chr; do
        samtools index -b "${OUT}/${pop}_${sample}_${chr}.bam"
    done < <(cut -f 2 ./Gene_locs.txt | sort -u)
}

export -f process_sample

while read -r s_line; do
    process_sample "$s_line" &
done < "$IN_NAME"
wait
