#!/bin/sh

TERASEQ="/root/TERA-Seq_manuscript"

for t in 5TERA 5TERA-short 5TERA3 TERA3 Akron5Seq mouse_SIRV dRNASeq RNASeq RiboSeq
do
	echo "downloading for $t"
	"$TERASEQ"/"$t"/download.sh
done
