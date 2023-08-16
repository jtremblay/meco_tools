#!/usr/bin/env Rscript

# Function that takes multiples gene count bed files and 
# merges them into on file (abundance matrix). Intended for the miRNA pipeline.
# Author: Julien Tremblay - jtremblay514@gmail.com

options(stringsAsFactors = FALSE)
mergeTablesMiRNA <- function(infiles, outfileRaw, outfileCpm, outfileRpkm, keepMatureOnly=TRUE) {

  library(edgeR)
  library(data.table)

  files = unlist(strsplit(infiles, ","))
  
  df = NULL
  
  tData = data.frame(fread(files[1], sep="\t", header=FALSE, showProgress=FALSE), check.names=FALSE)
  tData = tData[tData$V3 %in% c("miRNA", "miRNA_primary_transcript"),]
  print(head(tData))

  symbol_index = 9
  symbol = as.character(tData[,symbol_index])
  symbol_ref = gsub("^.*?;Name=(\\S+;*)Derives.*$", "\\1" , symbol)
  symbol_ref = gsub("^.*?;Name=(\\S+;*)$", "\\1" , symbol_ref)
  symbol_ref = gsub(";", "", symbol_ref)

  symbol_mirbase = gsub("^ID=(\\w+);.*?$", "\\1", symbol)
  symbol = NULL
  symbol = paste0(symbol_ref, "_", symbol_mirbase)
 
  df = cbind(df, symbol)

  #all lengths
  length_vec = abs(tData$V4 - tData$V5)

  #only mature mirna lengths
  tData_mature = tData[grepl("MIMAT", tData$V9),]
  length_vec_mature = abs(tData_mature$V4 - tData_mature$V5)
  colNames = c("feature_id")
  
  for(i in 1:length(files)){
    print(paste0("Reading file ", files[i]))
    currFile = files[i];
    currName = basename(currFile)
    currName = gsub(".cov", "", currName)
    colNames = append(colNames, currName)
    currTable = data.frame(fread(files[i], sep="\t", header=FALSE, showProgress=FALSE), check.names=FALSE)
    print(head(currTable))
    currTable = currTable[currTable$V3 %in% c( "miRNA", "miRNA_primary_transcript"),]
    rawCount = currTable[,ncol(currTable)]
 
    df = cbind(df, rawCount)
  }
  df = df[!duplicated(df), ]
  df = data.frame(df)
  names(df) = colNames
  df[is.na(df)] = 0
 
  if(keepMatureOnly == TRUE){
    df = df[grepl("MIMAT", df$feature_id),]
  }
  
  write.table(df, file=outfileRaw, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
 
  ## edgeR. Load table, remove low cpm. 
  df2 = df
  row.names(df2) = df2$feature_id
  df2$feature_id = NULL
  for(i in 1:ncol(df2)){
    df2[,i] = as.numeric(df2[,i])
  }

  # remove samples with only zeros, otherwise will crash
  df2 = df2[, apply(df2, 2, sum)!=0] 

  y <- DGEList(df2, remove.zeros=FALSE)
  y <- calcNormFactors(y, method="TMM")
  cpms = cpm(y)
  cpms = round(cpms, digits=3)
  cpms = as.data.frame(cpms)
  orig_colnames = colnames(cpms)
  cpms$feature_id = row.names(cpms)
  cpms = cpms[,c("feature_id", orig_colnames)]
  fwrite(cpms, outfileCpm, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
 
  y = NULL 
  y <- DGEList(df2, remove.zeros=FALSE)
  y <- calcNormFactors(y, method="TMM")
  rpkms = rpkm(y, gene.length=length_vec_mature)
  rpkms = round(rpkms, digits=3)
  rpkms = as.data.frame(rpkms)
  orig_colnames = colnames(rpkms)
  rpkms$feature_id = row.names(rpkms)
  rpkms = rpkms[,c("feature_id", orig_colnames)]
  fwrite(rpkms, outfileRpkm, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
  
  print("Done merging and normalizing tables...")
} 

usage=function(errM) {
  cat("\nUsage : Rscript mergeGeneAbundanceMiRNA.R [option] <Value>\n")
  cat("       -i        : List of files separated by a comma\n")
  cat("       -o        : outfile raw counts\n")
  cat("       -c        : outfile CPMs\n")
  cat("       -r        : outfile RPKMs\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 4) {
  usage("missing arguments")
}

skipNorm = NULL
## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-i") {
    infiles = ARG[i+1]
  } else if (ARG[i] == "-o") {
    outfile = ARG[i+1]
  } else if (ARG[i] == "-c") {
    outfileCpm = ARG[i+1]
  } else if (ARG[i] == "-r") {
    outfileRpkm = ARG[i+1]
  }
}

mergeTablesMiRNA(infiles, outfile, outfileCpm, outfileRpkm)
