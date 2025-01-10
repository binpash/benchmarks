# while read -r line
# do
#     cat $line |
#         iconv -c -t ascii//TRANSLIT |
#         pandoc +RTS -K64m -RTS --from html --to plain --quiet
# done

while read -r line; do
    iconv -c -t ascii//TRANSLIT < "$line" |
        pandoc +RTS -K64m -RTS --from html --to plain --quiet
done
