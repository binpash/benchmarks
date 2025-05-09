# #!/bin/bash 

# Calculate the bigrams (based on 4_3.sh script)
IN=${IN:-$SUITE_DIR/inputs/pg}
OUT=${1:-$SUITE_DIR/outputs/8.2_2/}
ENTRIES=${ENTRIES:-1000}
mkdir -p "$OUT"

pure_func() {
    input=$1
    IN=$2
    OUT=$3
    TEMPDIR=$(mktemp -d)

    cat "$IN/$input" | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$" > "${TEMPDIR}/${input}.input.words"
    tail +2 "${TEMPDIR}/${input}.input.words" > "${TEMPDIR}/${input}.input.nextwords"
    paste "${TEMPDIR}/${input}.input.words" "${TEMPDIR}/${input}.input.nextwords" | sort | uniq -c > "${TEMPDIR}/${input}.input.bigrams"
    awk '$1 == 2 {print $2, $3}' "${TEMPDIR}/${input}.input.bigrams" > "${OUT}/${input}.out"

    rm -rf "${TEMPDIR}"
}

export -f pure_func

ls "${IN}" | head -n "${ENTRIES}" | parallel -j "$(nproc)" pure_func {} "${IN}" "${OUT}"


echo "done"
