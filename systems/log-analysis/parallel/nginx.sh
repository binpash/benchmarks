#!/bin/bash
# Using GNU parallel:

mkdir -p $2

pure_func() {
    tempfile=$(mktemp)
    outputdir=$(mktemp -d)

    tee "$tempfile" | cut -d '"' -f3 | cut -d ' ' -f2 | sort | uniq -c | sort -rn > "$outputdir/output1.txt"

    parallel ::: \
        "awk '{print \$9}' $tempfile | sort | uniq -c | sort -rn > $outputdir/output2.txt" \
        "awk '(\$9 ~ /404/)' $tempfile | awk '{print \$7}' | sort | uniq -c | sort -rn > $outputdir/output3.txt" \
        "awk '(\$9 ~ /502/)' $tempfile | awk '{print \$7}' | sort | uniq -c | sort -r > $outputdir/output4.txt" \
        "awk -F'\"' '(\$2 ~ \"/wp-admin/install.php\"){print \$1}' $tempfile | awk '{print \$1}' | sort | uniq -c | sort -r > $outputdir/output5.txt" \
        "awk '(\$9 ~ /404/)' $tempfile | awk -F'\"' '(\$2 ~ \"^GET .*.php\")' | awk '{print \$7}' | sort | uniq -c | sort -r | head -n 20 > $outputdir/output6.txt" \
        "awk -F'\"' '{print \$2}' $tempfile  | awk '{print \$2}' | sort | uniq -c | sort -r > $outputdir/output7.txt" \
        "awk -F'\"' '(\$2 ~ \"ref\"){print \$2}' $tempfile | awk '{print \$2}' | sort | uniq -c | sort -r > $outputdir/output8.txt"

    # Combine results
    cat "$outputdir"/output*.txt
    rm -r "$outputdir" "$tempfile"
} 

# We parallelized the inner workload as it demonstrated a speedup.
# Since the number of files processed is fewer than the available cores
# inner parallelism resulted in improved performance.
# The inner parallelism is not necessary, but it is beneficial in this case.
# The commented out pure_func could also be used depending on the system.

export -f pure_func

export OUTPUT_DIR=$2

find "$1" -type f | parallel --jobs "$(nproc)" "cat {} | pure_func > $OUTPUT_DIR/{/}"