# count the packet number in a pcap file
tcpdump -nn -r ${INPUT2} | wc -l
