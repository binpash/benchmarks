#!/usr/bin/env Rscript
#
# A small script to subset two big tables
# Definitely faster than grep
#

library(optparse)
suppressMessages(library(dplyr))

# Read command line options and arguments
option_list <- list(
  make_option(
    c("-i", "--ifile"), type = "character",
    help = "Input table to subset from. Can be stdin (use: stdin). Quotes are ignored (quote=\"\").", metavar = "File"),
  make_option(
    c("-r", "--nonames"), action = "store_true", default = FALSE,
    help = "Input file (ifile) doesn't have column names."),  
  make_option(
    c("-e", "--exclude"), action = "store_true", default = FALSE,
    help = "Exclude the list instead of including it. Default is INCLUDE."),
  make_option(
    c("-c", "--column"),
    help = "Column name or column number to subset by (in ifile); must be column number if --nonames"),
  make_option(
    c("-l", "--ilist"), type = "character",
    help = "Input list to subset by (one per line, without header). Quotes are ignored (quote=\"\").", metavar = "File"),
  make_option(
    c("-o", "--ofile"), type = "character",
    help = "Output subset table. If empty output is stdout.", metavar = "File")
)
opt = parse_args(OptionParser(option_list = option_list))

### Testing variables
# print("USING TESTING VARIABLES!!!")
# opt<-NULL
# opt$ifile<-"/home/jan/projects/mourelatos11/projects/ribothrypsis/analysis/polya/results_example/20180620_1930_hsa-dRNASeq-HeLa-Total-REL3-1/rna_tails.tsv" # Table to subset from
# opt$ifile<-"test.sam"
# opt$column<-"read_id"
# opt$nonames<-T
# opt$column<-3
# opt$ilist<-"/home/jan/projects/mourelatos11/projects/ribothrypsis/analysis/polya/results_example/20180620_1930_hsa-dRNASeq-HeLa-Total-REL3-1/test.txt" # Table to subset by
# opt$ilist<-"/home/joppelt/projects/ribothrypsis/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.2/fastq/reads.1.sanitize.w_rel5.names.txt"
# opt$ofile<-"/home/joppelt/projects/mourelatos11/projects/ribothrypsis/analysis/polya/results_example/20180620_1930_hsa-dRNASeq-HeLa-Total-REL3-1/rna_tails.rel3.tsv" # Result table
# opt$ofile<-""
# print("USING TESTING VARIABLES!!!")
### Testing variables

#if(opt$ifile == "stdin"){
#  opt$ifile<-stdin()
#}

if(opt$nonames == TRUE){
  opt$column<-as.numeric(opt$column) # Such an easy source of errors...
  df_main<-read.delim(opt$ifile, header = FALSE, quote = "", sep="\t")  
}else{
  df_main<-read.delim(opt$ifile, quote = "", sep="\t")
}

df_subset<-read.table(opt$ilist, quote = "", sep="\t")

if(opt$exclude==TRUE){
        df_out <- df_main[!(df_main[, opt$column] %in% df_subset$V1), ]
}else{
	df_out <- df_main[df_main[, opt$column] %in% df_subset$V1, ]
}


if(opt$nonames == TRUE){
	if(length(opt$ofile)==0){
		write.table(x = df_out, file = stdout(), sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)
	}else{
		write.table(x = df_out, file = opt$ofile, sep="\t", row.names=FALSE, quote = FALSE, col.names = FALSE)
	}
}else{
	if(length(opt$ofile)==0){
		write.table(x = df_out, file = stdout(), sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)
	}else{
		write.table(x = df_out, file = opt$ofile, sep="\t", row.names=FALSE, quote = FALSE, col.names = TRUE)
	}
}
