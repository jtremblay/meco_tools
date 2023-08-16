#!/usr/bin/env Rscript

# Function that takes multiples gene count bed files and 
# merges them into on file (matrix). Intended for classic RNA-Seq
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca

options(stringsAsFactors = FALSE)
mergeHtseq <- function(infiles_htseq, outfileRaw, outfileNormalized) {
  library(edgeR)
  library(tools)
  #infiles_htseq = "~/R/ctrl-IC-1.cov,~/R/ctrl-IC-2.cov,~/R/ctrl-IC-3.cov"
  #outfileRaw = "~/R/outfile_raw.tsv"
  #outfileNormalized = "~/R/outfile_cpm.tsv"
  
  files_htseq = unlist(strsplit(infiles_htseq, ","))
  
  # symbols will be the same in both .htseq and .cov
  tData = read.table(files_htseq[1], sep="\t", header=FALSE, comment.char="")
  symbol = as.character(tData[,1])
  symbol = gsub("^id", "", symbol)

  df = NULL  
  df = cbind(df, symbol)
  df = data.frame(df)
  #df = df[grep("__", df[,1], invert=TRUE),]
  df = data.frame(df)
  df = df[order(df),]
  df = data.frame(df)
  
  colnames(df) = c("gene_id")
  
  for(i in 1:length(files_htseq)){
    print(paste0("Reading file ", files_htseq[i]))
    currFileHtseq = files_htseq[i];
    currName = file_path_sans_ext(basename(currFileHtseq))
    currName = gsub(".htseq", "", currName)
    currName = gsub(".cov", "", currName)
    tData = read.table(files_htseq[i], sep="\t", header=FALSE, comment.char="")
    colnames(tData) = c("gene_id", currName)
    tData$gene_id = gsub("^id", "", tData$gene_id)
    
    df = merge(df, tData, by="gene_id")
  }
  
  row.names(df) = df$gene_id
  df$gene_id = NULL
  
  y <- DGEList(counts=df, remove.zeros=FALSE)
  y_cpm = cpm(y)
  y_cpm = round(y_cpm, digits = 3)
  y_cpm = data.frame(y_cpm)
  
  #df[is.na(df)] = 0
  #dfNormalized[is.na(dfNormalized)] = 0
  colNames = colnames(df)
  df$gene_id = row.names(df)
  colNames = c("gene_id", colNames)
  df = df[ , colNames]
  
  colNames = colnames(y_cpm)
  y_cpm$gene_id = row.names(y_cpm)
  colNames = c("gene_id", colNames)
  y_cpm = y_cpm[ , colNames]
  
  
  write.table(df, file=outfileRaw, sep="\t", quote=FALSE, row.names=FALSE, col.name=TRUE)
  write.table(y_cpm, file=outfileNormalized, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  print("Done merging htseq tables...")
  
} 

usage=function(errM) {
  cat("\nUsage : Rscript mergeHtseqAbundance.R [option] <Value>\n")
  cat("       -i        : List of htseq files separated by a comma\n")
  cat("       -o        : outfile raw counts\n")
  cat("       -n        : outfile normalized\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 4) {
  usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-i") {
    infiles_htseq=ARG[i+1]
  } else if (ARG[i] == "-o") {
    outfile=ARG[i+1]
  } else if (ARG[i] == "-n") {
    outfileNormalized=ARG[i+1]
  }   
}

mergeHtseq(infiles_htseq, outfile, outfileNormalized)
