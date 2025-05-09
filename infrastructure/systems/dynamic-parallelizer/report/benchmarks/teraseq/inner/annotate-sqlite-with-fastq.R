#!/usr/bin/env Rscript
#
# Add (=annotate) sqlite database from fastq
#
# TODO: Take a multiple tables and an input instead of doing that one by one
#

library(RSQLite)
suppressMessages(library(dplyr))
library(optparse)

# Read command line options and arguments
option_list <- list(
  make_option(
    c("-d", "--database"), type = "character",
    help = "Input table database (sqlite) to annotate"),
  make_option(
    c("-n", "--db_col_bind"), type = "character",
    help = "Database name column to annotate by (~overlap with the list)"),
  make_option(
    c("-c", "--db_col_add"), type = "character",
    help = "Name of the column to add to the database"),
  make_option(
    c("-t", "--db_tables"), type = "character",
    help = "Name of ONE (so far) database table to annotate"),
  make_option(
    c("-i", "--ifile"), type = "character",
    help = "Name of the input fastq to annotate by")
)
opt = parse_args(OptionParser(option_list = option_list))

### Testing variables
#opt<-NULL
#opt$db_col_add <- "rel3"
#opt$db_col_bind <- "qname"
#opt$db_tables <- "genome" # [genome, transcr, ribo]
#opt$database<-"sqlite.db"
#opt$ifile<-"reads.1.sanitize.w_rel3.fastq.gz" # input fastq
### Testing variables

if(length(grep(".gz$", opt$ifile))==1){
  annot_list <- read.table(pipe(paste("zcat", opt$ifile, "| paste - - - - | cut -f1 | sed 's/^@//g'")), header=F)
}else{
  annot_list <- read.table(pipe(paste("cat", opt$ifile, "| paste - - - - | cut -f1 | sed 's/^@//g'")), header=F)
}

dbcon <- dbConnect(RSQLite::SQLite(), opt$database)

# https://stackoverflow.com/questions/43967965/adding-column-to-sqlite-database
for(db_table in unique(opt$db_tables)){
  
  print(paste0("Annotating table: ", db_table, "."))
  
  db_match <- as.integer(as.logical(dbGetQuery(dbcon, paste("SELECT", opt$db_col_bind, "FROM", db_table))[, 1] %in% annot_list$V1)) # Get the overlap between the annotation table and the list to add to the annotation
  db_match[db_match==0]<-NA # Convert zeros to NA to fit to Manolis scheme (TRUE = 1, FALSE = NA)
  
  dbExecute(dbcon, paste("ALTER TABLE", db_table, "ADD COLUMN", opt$db_col_add, "INT(1)")) # Initialize new column

  # Prepare a table to make a merge (annotate) the sqlite  
  johndoe <- dbGetQuery(dbcon, paste("SELECT", opt$db_col_add, "FROM", db_table)) # Get the whole table - we need it for the "merge" afterwards
  johndoe <- data.frame(placeholder=db_match,
             id=rownames(johndoe))
  colnames(johndoe)[1]<-opt$db_col_add # Rename the first column to fit to the variable - probably easiest solution

  dbBegin(dbcon)
  dbExecute(dbcon, paste0("UPDATE ", db_table, " SET ", opt$db_col_add, "= :", opt$db_col_add, " where rowid = :id"),
            params=johndoe)
  dbCommit(dbcon)
}

# Check if we actually added something to the column) and that it fits to the input annotation list - we should see a match between the two rows
print("Check the annotation fits - we should see a match between the input annotation and the sqltable")
tmp.check <- head(dbGetQuery(dbcon, paste0("SELECT ", opt$db_col_add, ",", opt$db_col_bind, " FROM ", db_table)), n = 20)[, 1] # Convert back from 1,NA to 1,0
tmp.check[is.na(tmp.check)]<-0
tmp.check == as.integer(as.logical((head(dbGetQuery(dbcon, paste0("SELECT ", opt$db_col_add, ",", opt$db_col_bind, " FROM ", db_table)), n = 20)[, 2] %in% annot_list[, 1])))

dbDisconnect(dbcon) # Exit database and save