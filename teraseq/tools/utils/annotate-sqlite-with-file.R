#!/usr/bin/env Rscript
#
# Add (=annotate) sqlite database from file, optionally with int values
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
    c("-r", "--round"), action="store_true", default=FALSE,
    help = "Round up (to integer) values in the --ifile"),
  make_option(
    c("-i", "--ifile"), type = "character",
    help = "Name of the input file to annotate by")
)
opt = parse_args(OptionParser(option_list = option_list))

### Testing variables
#opt<-NULL
#opt$db_col_add <- "polya"
#opt$db_col_bind <- "qname"
#opt$db_tables <- "genome" # [genome, transcr, ribo]
#opt$database<-"db/sqlite.db.bckp"
#opt$ifile<-"align/reads.1.sanitize.noribo.toTranscriptome.sorted.polya.filt.txt" # input fastq
#opt$round<-TRUE
### Testing variables

if(length(grep(".gz$", opt$ifile))==1){
  annot_list <- read.table(pipe(paste("zcat", opt$ifile)), header=F, stringsAsFactors=F)
}else{
  annot_list <- read.table(opt$ifile, header=F, stringsAsFactors=F)
}

colnames(annot_list)[1]<-"read_id"

if(ncol(annot_list)==2){
  if(opt$round==TRUE){
    if(typeof(annot_list$V2)=="double"){
      annot_list$V2<-as.integer(round(annot_list$V2, 0))
    }else{
      print("Don't know how to round, the number is not number or is not double")
    }
  }
}

dbcon <- dbConnect(RSQLite::SQLite(), opt$database)

# https://stackoverflow.com/questions/43967965/adding-column-to-sqlite-database
for(db_table in unique(opt$db_tables)){

  print(paste0("Annotating table: ", db_table, "."))

  db_match<-as.data.frame(dbGetQuery(dbcon, paste("SELECT", "rowid,", opt$db_col_bind, "FROM", db_table))[, 1:2], stringsAsFactors=F)
  colnames(db_match)[2]<-"read_id"

  db_match$sortid<-1:nrow(db_match)

  if(ncol(annot_list)==2){
    db_match<-merge(db_match, annot_list, all.x=T)
#    db_match[is.na(db_match[,2]),2]<--1 # set NA values to -1
  }else{
    db_match$V3 <- as.integer(as.logical(db_match[, 2] %in% annot_list[, 1])) # Get the overlap between the annotation table and the list to add to the annotation
    db_match$V3[db_match$V3==0]<-NA # Convert zeros to NA to fit to Manolis scheme (TRUE = 1, FALSE = NA)
  }

  db_match<-db_match[order(db_match$sortid), ]
  db_match$sortid<-NULL

  dbExecute(dbcon, paste("ALTER TABLE", db_table, "ADD COLUMN", opt$db_col_add, "INT(10)")) # Initialize new column

  # Prepare a table to make a merge (annotate) the sqlite
  johndoe <- dbGetQuery(dbcon, paste("SELECT", opt$db_col_add, "FROM", db_table)) # Get the whole table - we need it for the "merge" afterwards
  johndoe <- data.frame(placeholder=db_match,
                          id=rownames(johndoe))

  colnames(johndoe)[3]<-opt$db_col_add # Rename the first column to fit to the variable - probably easiest solution
  johndoe$id<-johndoe$placeholder.id
  johndoe$placeholder.id<-NULL
  johndoe$placeholder.read_id<-NULL

  dbExecute(dbcon, paste0("UPDATE ", db_table, " SET ", opt$db_col_add, " = :", opt$db_col_add, " where rowid= :id"),
              params=johndoe)
}

# Check if we actually added something to the column) and that it fits to the input annotation list - we should see a match between the two rows
print("Check the annotation fits - we should see a match between the input annotation and the sqltable. Result should be a number and 0")
check <- dbGetQuery(dbcon, paste0("SELECT ", opt$db_col_bind, ",", opt$db_col_add, " FROM ", db_table))
if(ncol(annot_list==2)){
  colnames(annot_list)[1]<-opt$db_col_bind
  check<-merge(check, annot_list, all.x=T)
  print(sum(check[,2]==check[,3], na.rm=T))
  print(sum(check[,2]!=check[,3], na.rm=T))
}else{
  check<-check[!is.na(check[, opt$db_col_add]),]
  sum(annot_list$read_id %in% check[!is.na(check[, opt$db_col_add]), opt$db_col_bind])
  sum(annot_list$read_id %in% check[is.na(check[, opt$db_col_add]), opt$db_col_bind])
}

dbDisconnect(dbcon) # Exit database and save
