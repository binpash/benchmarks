#!/usr/bin/env Rscript
#
# Plot heatmap of real matches of adapter in libraries with various cutadapt settings
# Also checks if the random matches are 0
#
# Script has three arguments:
#	Directory with cutadapt logfiles
#	Output plot file (png)
#	tsv output of cutadapt-transcripts.R (random matches with transcripts)
#

args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library("dplyr"))
library("ggplot2")
# library("reshape2")
# library("stringr")

indir <- args[1] # "results/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/cutadapt/logfiles"
ofile <- args[2] # "results/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/cutadapt/cutadapt-heatmap.png"
trans_file <- args[3] # "results/transcripts/REL5.long/cutadapt-heatmap.tsv" # file with results of cutadapt-transcripts.R - random matches with transcripts

input_files <- list.files(path = indir, pattern = "cutadapt..*\\.log$", full.names = T)

origs<-grep(pattern = ".orig.", x = input_files, fixed = T)
shufs<-input_files[-origs]
origs<-input_files[origs]

full_tab <- data.frame(
  overlap = numeric(length(origs)),
  error = numeric(length(origs)),
  trimmed = numeric(length(origs))
  )

for (counter in 1:length(origs)) {
  input_file<-origs[counter]

  settings<-paste(strsplit(basename(input_file), ".", fixed = T)[[1]][3:4], collapse = ".")

  len <- strsplit(x = basename(input_file), split = ".", fixed = T)[[1]][3] %>% gsub(pattern = "l", replacement = "")
  err <- strsplit(x = basename(input_file), split = ".", fixed = T)[[1]][4] %>% gsub(pattern = "e", replacement = "")

  # Check if all the shuffled adapters are 0 and move to next it not
  shuf<-shufs[grep(pattern = settings, shufs, fixed = T)]

  trims<-0
  for(shuf_cur in shuf){
    shuf_trim <- system(paste0("grep \"Reads with adapters:\" ", shuf_cur, "| tr -s ' ' | cut -d ' ' -f 5"), intern = TRUE) # read only number of trimmed adapts
    shuf_trim <- as.numeric(stringr::str_replace_all(string = shuf_trim, pattern = "[()%]", replacement = ""))
    trims <- trims + shuf_trim
  }

  if(trims != 0){
    trimmed <- NA
  }else{
    trimmed <- system(paste0("grep \"Reads with adapters:\" ", input_file, "| tr -s ' ' | cut -d ' ' -f 5"), intern = TRUE) # read only number of trimmed adapts
    trimmed <- as.numeric(stringr::str_replace_all(string = trimmed, pattern = "[()%]", replacement = ""))
  }

  full_tab$overlap[counter] <- len
  full_tab$error[counter] <- err
  full_tab$trimmed[counter] <- trimmed
}

# https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/

# Ranges are open on the left instead of right (default)
m <- full_tab %>%
  mutate(countfactor = cut(trimmed,
                           breaks = c(-1, 0.001, 5, 10, 15, 20, 25, max(full_tab$trimmed, na.rm = T)),
                           labels = c("0", "0-5", "5-10", "10-15", "15-20", "20-25", ">25"),
						   right=F, include.lowest=T
                            )
        )

# Add text inside the tiles https://stackoverflow.com/questions/62007556/using-ggplot2-to-print-text-with-n-in-them
p <- ggplot() +
  geom_tile(data = m, aes(x = error, y = overlap, fill = countfactor), colour = "white", size = 0.2) +
  labs(x = "Maximum error rate", y = "Minimum overlap length") +
  scale_y_discrete(expand = c(0.015, 0)) +
  scale_fill_manual(values = c("#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da"), na.value = "grey90") +
  theme_classic(base_size = 8) +
  guides(fill = guide_legend(title = "Percentage of trimmed\nreads"))

if(file.exists(trans_file)){
  trans <- read.table(trans_file, sep="\t", row.names = 1)
  colnames(trans)<-trans[1 ,]
  trans <- trans[-1 ,]
  trans<-reshape2::melt(as.matrix(trans))
  colnames(trans) <- c("overlap","error","trimmed")
  trans <- trans %>%
    arrange(overlap, error)
  trans$overlap<-as.character(trans$overlap)
  trans$error<-as.character(trans$error)

  # add strikethrough text and bold and color
  trans$font <- "plain"
  trans$font[trans$trimmed<=0.2] <- "bold"
  trans$col <- "black"
  trans$col[trans$trimmed>0.2]<-"darkgrey"

  p <- p + geom_text(data = trans, aes(x = error, y = overlap, label = trimmed), fontface = trans$font, color =  trans$col, size=1.8)
}

ggsave(p, filename = ofile, height = 5.5, width = 8.8, units = "in", dpi = 200)

# Make a table for export
full_tab <- reshape2::dcast(data = m, formula = overlap ~ error, value.var = "trimmed") # unmelt
write.table(x = full_tab, file = paste0(sub(".[^.]+$", "", ofile), ".tsv"), row.names = F, sep = "\t", quote = F)
