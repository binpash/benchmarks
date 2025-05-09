#!/bin/sh
try=/srv/hs/deps/try/try

# Step 1: Quality Control with FastQC
$try -y echo "Running FastQC..."
mkdir -p fastqc_output
for sample in SRR10045016_1.fastq SRR10045017_1.fastq SRR10045018_1.fastq SRR10045019_1.fastq SRR10045020_1.fastq SRR10045021_1.fastq SRR10045016_2.fastq SRR10045017_2.fastq SRR10045018_2.fastq SRR10045019_2.fastq SRR10045020_2.fastq SRR10045021_2.fastq; do
    $$try -y ./FastQC/fastqc -t 2 data_chr19/${sample} -o fastqc_output/
done

# Step 2: Trimming with Trimmomatic
$try -y echo "Trimming sequences with Trimmomatic..."
$try -y mkdir -p trim_output
for sample in SRR10045016 SRR10045017 SRR10045018 SRR10045019 SRR10045020 SRR10045021; do
	    $try -y java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 2 data_chr19/${sample}_1.fastq data_chr19/${sample}_2.fastq \
		        -baseout trim_output/${sample}/ ILLUMINACLIP:trimmomaticAdapters/TruSeq3-PE-2.fa:2:30:10 TRAILING:30
done

# Step 3: Post-trimming Quality Control
$try -y echo "Running FastQC on trimmed sequences..."
mkdir -p fastqc_trimmed_output
for sample in SRR10045016_1P SRR10045017_1P SRR10045018_1P SRR10045019_1P SRR10045020_1P SRR10045021_1P SRR10045016_2P SRR10045017_2P SRR10045018_2P SRR10045019_2P SRR10045020_2P SRR10045021_2P; do
    $try -y ./FastQC/fastqc -t 2 trim_output/${sample} -o fastqc_trimmed_output/
done

# Step 4: Alignment with STAR
$try -y echo "Generating STAR index..."
$try -y mkdir -p STAR_index_chr19
$try -y ../STAR/source/STAR --runThreadN 2 --runMode genomeGenerate --genomeDir STAR_index_chr19/ --genomeFastaFiles genome/chr19.fa --sjdbGTFfile genome/chr19_Homo_sapiens.GRCh38.95.gtf --genomeSAindexNbases 11

echo "Aligning sequences with STAR..."
mkdir -p STAR_output
for sample in SRR10045016 SRR10045017 SRR10045018 SRR10045019 SRR10045020 SRR10045021; do
	$try -y mkdir STAR_output/${sample}
	$try -y ../STAR/source/STAR --runThreadN 2 --genomeDir STAR_index_chr19/ --readFilesIn trim_output/${sample}_1P trim_output/${sample}_2P --sjdbGTFfile genome/chr19_Homo_sapiens.GRCh38.95.gtf --outFileNamePrefix STAR_output/${sample}/
done

# Step 5: Counting with htseq-count
echo "Counting alignments with htseq-count..."
$try -y mkdir -p htseq_output
$try -y htseq-count -n 1 $(for sample in SRR10045016 SRR10045017 SRR10045018 SRR10045019 SRR10045020 SRR10045021; do try -y echo STAR_output/${sample}/Aligned.out.sam; done) genome/chr19_Homo_sapiens.GRCh38.95.gtf > htseq_output/SRR10045016-17-18-19-20-21_counts.csv

$try -y echo "RNA-Seq analysis completed."
