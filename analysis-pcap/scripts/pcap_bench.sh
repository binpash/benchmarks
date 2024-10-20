time tcpdump -nn -r <(cat ${INPUT}) -A 'port 53'| sort | uniq |grep -Ev '(com|net|org|gov|mil|arpa)'     

# Throws permission denied
# time tcpdump -nn -r ${INPUT} -A 'port 53'| sort | uniq |grep -Ev '(com|net|org|gov|mil|arpa)'     

time tcpdump -nn -r ${INPUT2} -A -c 1000000 | sort |uniq |grep -Ev '(com|net|org|gov|mil|arpa)'
