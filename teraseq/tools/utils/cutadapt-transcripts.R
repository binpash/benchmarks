#!/usr/bin/env Rscript
#
# Plot heatmap of random matches of adapter in transcripts with various cutadapt settings
#
# Script has two arguments:
#	Directory with cutadapt logfiles
#	Output plot file (png)
#

args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library("dplyr"))
library("ggplot2")
# library("reshape2")
# library("stringr")

indir <- args[1] # "results/transcripts/REL5.long/logfiles"
ofile <- args[2] # "results/transcripts/REL5.long/cutadapt-heatmap.png"

input_files <- list.files(path = indir, pattern = "cutadapt..*\\.log$", full.names = T)

full_tab <- data.frame(
  overlap = numeric(length(input_files)),
  error = numeric(length(input_files)),
  trimmed = numeric(length(input_files))
)

for (counter in 1:length(input_files)) {
  input_tab <- NULL
  input_file <- input_files[counter]

  len <- strsplit(x = basename(input_file), split = ".", fixed = T)[[1]][2] %>% gsub(pattern = "l", replacement = "")
  err <- strsplit(x = basename(input_file), split = ".", fixed = T)[[1]][3] %>% gsub(pattern = "e", replacement = "")

  # Check if the line we want exists
  start_line <- system(paste("grep -nr \"Reads with adapters: \"", input_file, "| cut -d \":\" -f1"), intern = TRUE) # Get line with no of trimmed sequences
  start_line <- as.numeric(start_line)

  if (length(start_line) != 0) { # Skip empty log files
    trimmed <- system(paste0("grep \"Reads with adapters:\" ", input_file, "| tr -s ' ' | cut -d ' ' -f 5"), intern = TRUE) # read only number of trimmed adapts
    trimmed <- as.numeric(stringr::str_replace_all(string = trimmed, pattern = "[()%]", replacement = ""))
  } else {
    trimmed <- NA
  }

  full_tab$overlap[counter] <- len
  full_tab$error[counter] <- err
  full_tab$trimmed[counter] <- trimmed
}

# https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/

# Ranges are open on the left instead of right (default)
m <- full_tab %>%
  mutate(countfactor = cut(trimmed,
    breaks = c(-1, 0.2, 1, 2, 5, 10, 25, max(full_tab$trimmed, na.rm = T)),
    labels = c("<0.2", "0.2-1", "1-2", "2-5", "5-10", "10-25", ">=25"),
	right=F, include.lowest=T
  )) %>%
  mutate(countfactor = factor(as.character(countfactor), levels = rev(levels(countfactor))))

p <- ggplot(m, aes(x = error, y = overlap, fill = countfactor)) +
  geom_tile(colour = "white", size = 0.2) +
  labs(x = "Maximum error rate", y = "Minimum overlap length") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c("#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da"), na.value = "grey90") +
  theme_classic(base_size = 8) +
  guides(fill = guide_legend(title = "Percentage of trimmed\ntranscripts"))

ggsave(p, filename = ofile, height = 5.5, width = 8.8, units = "in", dpi = 200)

# Make a table for export
full_tab <- reshape2::dcast(data = m, formula = overlap ~ error, value.var = "trimmed") # unmelt
write.table(x = full_tab, file = paste0(sub(".[^.]+$", "", ofile), ".tsv"), row.names = F, sep = "\t", quote = F)
