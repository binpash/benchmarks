#!/bin/sh
# showfile--Shows the contents of a file, including additional useful info
IFS='
'

width=72
for input in "$@"
do
    lines=$(wc -l < "$input" | tr -d ' ')
    chars=$(wc -c < "$input" | tr -d ' ')
    owner=$(ls -ld "$input" | awk '{print $3}')
    echo "-----------------------------------------------------------------"
    echo "File $input ($lines lines, $chars characters, owned by $owner):"
    echo "-----------------------------------------------------------------"
    
    file_content=$(cat "$input")
    for line in $(echo "$file_content"); do
        if [ ${#line} -gt $width ]; then
            formatted_line=$(echo "$line" | fmt)
            echo "$formatted_line" | {
                first_line=true
                for subline in $(cat); do
                    if $first_line; then
                        echo " $subline"
                        first_line=false
                    else
                        echo "+ $subline"
                    fi
                done
            }
        else
            echo " $line"
        fi
    done
    echo "-----------------------------------------------------------------"
done
