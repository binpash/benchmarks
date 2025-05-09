#!/bin/bash

RESOURCE_DIR="$RESOURCE_DIR"
INPUT_DIR="$INPUT_DIR"
MAPPING_FILE="$MAPPING_FILE"
TOOL="$TOOL"
OUTPUT_DIR="$OUTPUT_DIR"


# Process each compressed log file in the input directory
for gz_file in $(find "${RESOURCE_DIR}" -type f -name '*.gz'); do

    # Determine the base name without the .gz extension
    INPUT=$(basename "${gz_file%.gz}")
    input_file="$OUTPUT_DIR/$INPUT.log"
    response_codes_file="$OUTPUT_DIR/${INPUT}_response_codes.log"
    broken_links_file="$OUTPUT_DIR/${INPUT}_404_broken_links.log"
    errors_file="$OUTPUT_DIR/${INPUT}_502_errors.log"
    requests_wp_admin_file="$OUTPUT_DIR/${INPUT}_requests_wp_admin.log"
    php_files_file="$OUTPUT_DIR/${INPUT}_404_php_files.log"
    most_requested_urls_file="$OUTPUT_DIR/${INPUT}_most_requested_urls.log"
    most_requested_urls_ref_file="$OUTPUT_DIR/${INPUT}_most_requested_urls_ref.log"

    # Decompress and convert the binary log to ASCII log
    gzip -dc "$gz_file" | "$TOOL" "$MAPPING_FILE" > "$input_file"

    # Sort by response codes
    cat "$input_file" | cut -d "\"" -f3 | cut -d ' ' -f2 | sort | uniq -c | sort -rn > "$response_codes_file"

    # Find broken links (404)
    awk "(\$9 ~ /404/)" "$input_file" | awk "{print \$7}" | sort | uniq -c | sort -rn > "$broken_links_file"

    # For 502 (bad-gateway)
    awk "(\$9 ~ /502/)" "$input_file" | awk "{print \$7}" | sort | uniq -c | sort -r > "$errors_file"

    # Who are requesting broken links (or URLs resulting in 502)
    awk -F\" "(\$2 ~ \"/wp-admin/install.php\"){print \$1}" "$input_file" | awk "{print \$1}" | sort | uniq -c | sort -r > "$requests_wp_admin_file"

    # 404 for php files - mostly hacking attempts
    awk "(\$9 ~ /404/)" "$input_file" | awk -F\" "(\$2 ~ \"^GET .*\.php\")" | awk "{print \$7}" | sort | uniq -c | sort -r | head -n 20 > "$php_files_file"

    # Most requested URLs
    awk -F\" "{print \$2}" "$input_file" | awk "{print \$2}" | sort | uniq -c | sort -r > "$most_requested_urls_file"

    # Most requested URLs containing "ref"
    awk -F\" "(\$2 ~ \"ref\"){print \$2}" "$input_file" | awk "{print \$2}" | sort | uniq -c | sort -r > "$most_requested_urls_ref_file"

done
