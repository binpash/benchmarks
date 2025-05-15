#!/usr/bin/env Rscript
#
# Fix/add transcript field based on exons (start, end, strand)
# Checks if transcript feature exists, compares it with exon-based info and replaces it if it's different
# Doesn't consider UTRs right now
#
# Original made for SIRV gtf annotation which has often incorrect strand in transcript field
#

# TODO: add UTRs and start/stop codons to transcript start/end coords
# TODO: if we don't have all transcript lines check those which exist and add missing

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))

args <- commandArgs(trailingOnly=TRUE)

ingtf<-args[1]
outgtf<-args[2]

#ingtf<-"E:/Results/zissimos/spikeins/sirv/SIRVome_isoforms_C_170612a.E2.gtf"
#outgtf<-"E:/Results/zissimos/spikeins/sirv/SIRVome_isoforms_C_170612a.E2.fixed.gtf"

gtf<-read.table(ingtf, sep="\t", header=F, stringsAsFactors = F, quote = "")

gtf$transcript_id<-lapply(gtf$V9, function(x) stringr::str_split(x, "; | ")[[1]][4]) %>%
  unlist() %>%
  gsub(pattern = ";", replacement = "") %>%
  gsub(pattern = "\"", replacement = "") # transcript_id

trans<-gtf[gtf$V3=="transcript",]
exons<-gtf[gtf$V3=="exon",]

# Get the last exon ONLY for + transcripts
starts<-exons %>%
  group_by(transcript_id) %>%
  summarise(start=min(V4))

# Get the last exon ONLY for - transcripts
ends<-exons %>%
  group_by(transcript_id) %>%
  summarise(end=max(V5))

# Get strands
strands<-unique(exons[,c(7, ncol(exons))])

trans_mod<-starts %>%
  left_join(ends) %>%
  left_join(strands) %>%
  rename(strand="V7")

if(length(unique(gtf$transcript_id))!=nrow(trans)){
  print("Whoops, number of transcript feature lines is not the same as number of transcripts in the gtf, removing all and making new ones")

  exons$V9<-lapply(exons$V9, function(x) stringr::str_split(x, " exon_assignment")[[1]][1]) %>% unlist()

  trans_mod<-exons %>%
    select(V1, V2, V3, V6, V8, V9, transcript_id) %>%
    unique() %>%
    left_join(trans_mod) %>%
    mutate(V3="transcript") %>%
    rename(V4="start", V5="end", V7="strand")

}else{
  print("Number of transcript feature lines is the same as number of transcripts in the gtf, checking existing, comparing with exon-extracted info and replacing if they are not the same.")

  # Merge exon-extracted info with the original trans
  trans_mod<-trans_mod %>%
    right_join(trans) %>%
    as.data.frame()

  # Get trans lines which are different
  trans_mod[paste(trans_mod[, c("start", "end", "strand")], collapse = ".")==paste(trans_mod[, c("V4", "V5", "V7")], collapse = ".")]

  # Add helper columns for comparison
  trans_mod <- trans_mod %>%
    group_by(transcript_id) %>%
    mutate(tag1 = paste(c(V4, V5, V7), collapse=".")) %>%
    mutate(tag2 = paste(c(start, end, strand), collapse="."))

  # Get lines where the transcript info is different from the one from exons
  print("These transcript lines are different between the original annotation and exon-extracted info. Please check what's wrong.")
  print("Vx named columns are from the original gtf while the 'normal' named are extracted from exons.")
  index<-trans_mod$tag1!=trans_mod$tag2
  trans_mod[index,] %>%
    select(-tag1, -tag2) %>%
    print()
  #trans_mod[trans_mod$tag1==trans_mod$tag2,] %>% select(-tag1, -tag2) # The same from transcript and exons

  # Replace wrong transcript lines
  trans_mod[index, c("V4", "V5", "V7")]<-trans_mod[index, c("start", "end", "strand")]
}

gtf_out <- trans_mod %>%
  ungroup() %>%
  select(V1, V2, V3, V4, V5, V6, V7, V8, V9) %>%
  rbind(gtf %>% select(-transcript_id) %>% filter(V3!="transcript")) %>%
  arrange(V9, V4, V5)

write.table(gtf_out, file = outgtf, sep="\t", quote = F, row.names=F, col.names = F)
