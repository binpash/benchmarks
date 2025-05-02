#!/bin/bash

_5TERA="\
hsa.dRNASeq.HeLa.polyA.REL5OH.long.1/fastq/reads.1.fastq.gz \
hsa.dRNASeq.HeLa.polyA.REL5.long.1/fastq/reads.1.fastq.gz \
hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/fastq/reads.1.fastq.gz \
hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1/fastq/reads.1.fastq.gz
"
_5TERA_short="\
hsa.dRNASeq.HeLa.polyA.REL5.1/fastq/reads.1.fastq.gz \
hsa.dRNASeq.HeLa.polyA.PNK.REL5.1/fastq/reads.1.fastq.gz
"

_5TERA3="\
hsa.dRNASeq.HeLa.total.REL5.long.REL3.4/fastq/reads.1.fastq.gz \
hsa.dRNASeq.HeLa.total.REL5.long.REL3.5/fastq/reads.1.fastq.gz \
hsa.dRNASeq.HeLa.total.REL5.long.REL3.6/fastq/reads.1.fastq.gz
"

_TERA3="\
hsa.dRNASeq.HeLa.total.REL3.1/fastq/reads.1.fastq.gz \
hsa.dRNASeq.HeLa.total.REL3.2/fastq/reads.1.fastq.gz \
hsa.dRNASeq.HeLa.total.REL3.3/fastq/reads.1.fastq.gz
"

_Akron5Seq="\
hsa.Akron5Seq.HeLa.whole.2/fastq/reads.1.fastq.gz \
hsa.Akron5Seq.HeLa.whole.2/fastq/reads.2.fastq.gz
"

_mouse_SIRV="\
mmu.dRNASeq.inclSIRV.PRJEB27590.ERR3363659.1/fastq/reads.1.fastq.gz \
mmu.dRNASeq.inclSIRV.PRJEB27590.ERR2680379.1/fastq/reads.1.fastq.gz \
mmu.dRNASeq.inclSIRV.PRJEB27590.ERR3363657.1/fastq/reads.1.fastq.gz \
mmu.dRNASeq.inclSIRV.PRJEB27590.ERR2680375.1/fastq/reads.1.fastq.gz
"

_dRNASeq="hsa.dRNASeq.HeLa.polyA.1/fastq/reads.1.fastq.gz"

_RiboSeq="hsa.RiboSeq.HeLa.async.2/fastq/reads.1.fastq.gz"

_RNASeq="\
hsa.RNASeq.HeLa.xxx.polyA.ENCSR000CPR.1/fastq/reads.1.fastq.gz \
hsa.RNASeq.HeLa.xxx.polyA.ENCSR000CPR.1/fastq/reads.2.fastq.gz
"

for f in $_5TERA $_5TERA_short $_5TERA3 $_TERA3 $_Akron5Seq $_mouse_SIRV $_dRNASeq $_RiboSeq $_RNASeq
do
	url="https://atlas.cs.brown.edu/data/teraseq/$f"
	path="samples/$f"
	
	echo "downloading $path"
	echo "url: $url"

	mkdir -p $(dirname $path)
	if [[ ! -f $path ]]
	then
		wget --no-check-certificate $url -O $path
	else
		echo "already exists. skipping..."
	fi
done
