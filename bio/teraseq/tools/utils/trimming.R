#!/usr/bin/env Rscript
#
# Plot adapter trimming
# Script has two arguments:
#	Input TSV file with number of trimmed reads in the samples with length\tsamp1_trim\tsamp2_trim\t...
#	Output PDF
#

args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggplot2"))

norm1 <- function(x){(x-min(x))/(max(x)-min(x))} # normalized to one function

trimming<-read.table(args[1], sep="\t", header=T)

colnames(trimming)[colnames(trimming)=="length"]<-"Length"
colnames(trimming)[colnames(trimming)=="CIP.decap"]<-"Cap-Poly(A)"
colnames(trimming)[colnames(trimming)=="decap"]<-"Cap & 5P-Poly(A)"
colnames(trimming)[colnames(trimming)=="X5P"]<-"5P-Poly(A)"
colnames(trimming)[colnames(trimming)=="X5OH"]<-"5OH-Poly(A)"

trimming.m<-reshape2::melt(trimming, id.vars="Length")

trimming.m<-trimming.m %>%
  group_by(variable) %>%
  mutate(freq = value / sum(value)) %>%
  rename("Library"="variable")

max_val<-trimming.m %>%
  group_by(Library) %>%
  top_n(1, freq) %>%
  select(Length, Library)

breaks = c(seq(20, 80, by=10), median(max_val$Length), 58) # generate break positions
labels = as.character(breaks) # and labels

p <- ggplot(trimming.m, aes(x=Length, y=freq, color=Library)) +
  geom_line() +
  scale_x_continuous(limits = c(20, 80), breaks = breaks, labels = labels) +
  theme_classic() +
  theme(legend.position="bottom") +
#  ggtitle("Removed sequence length") +
  xlab("Removed length") + ylab("Frequency")

p <- p +
#  geom_vline(data=max_val, aes(xintercept=as.numeric(Length), color=Library)) +
  geom_vline(xintercept=median(max_val$Length), color="orange") +
  geom_vline(xintercept = 58, color="grey") +
  theme(plot.title = element_text(hjust = 0.5))

pdf(args[2])
  print(p)
dev.off()
