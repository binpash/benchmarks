#!/bin/bash

tests=(
	"5TERA"
	"5TERA-short"
	"5TERA3"
	"TERA3"
	"dRNASeq"
	"RNASeq"
	"RiboSeq"
	"Akron5Seq"
	"mouse_SIRV"
)

out="$PWD/mem_usage.txt"

for t in ${tests[@]}
do
	cd $t
	pwd
	date | tee -a $out
	echo "testing $t" | tee -a $out
	container_id=$(docker run -d --rm --privileged --name teraseq_mem_usage -v .:/root/TERA-Seq_manuscript/samples teraseq20-data-hs /bin/bash /root/TERA-Seq_manuscript/samples/run.sh)
	echo "started printer"
	while true
	do
		if docker ps -f id=$container_id | grep teraseq_mem_usage >/dev/null
		then
			date >> $out
			docker stats --no-stream teraseq_mem_usage 2>/dev/null | grep teraseq_mem_usage >> $out
			sleep 15
		else
			break
		fi
	done
	echo "stopped printer"
	cd ..
done

