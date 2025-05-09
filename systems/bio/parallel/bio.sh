#!/bin/bash


cd "$(dirname "$0")"/.. || exit 1

IN="$1"
IN_NAME="$2"
OUT="$3"

process_sample() {
    sample=$(echo "$s_line" | cut -d " " -f 2)
    pop=$(echo "$s_line" | cut -f 1 -d " ")
    link=$(echo "$s_line" | cut -f 3 -d " ")

    echo "Processing Sample $sample"

    samtools view -H "${IN}/$sample".bam | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' \
      | samtools reheader - "${IN}/$sample".bam > "${OUT}/$sample"_corrected.bam

    samtools index -b "${out_dir}/$sample"_corrected.bam

    cut -f 2 ./Gene_locs.txt | sort | uniq | while read chr; do
      echo "Isolating Chromosome $chr from sample ${OUT}/$sample"
      samtools view -b "${OUT}/$sample"_corrected.bam chr"$chr" > "${OUT}/$pop"_"$sample"_"$chr".bam
      echo "Indexing Sample ${OUT}/$pop_$sample_$chr"
      samtools index -b "${OUT}/$pop"_"$sample"_"$chr".bam
    done;
    # Optionally remove intermediate files
    # rm "${out_dir}/$sample"_corrected.bam
    # rm "${out_dir}/$sample"_corrected.bam.bai
}
export -f process_sample

export IN OUT

# Parallelize the processing of samples
cat "$IN_NAME" | parallel --jobs "$(nproc)" process_sample {} "$IN" "$OUT"