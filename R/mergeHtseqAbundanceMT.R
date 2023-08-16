#!/usr/bin/env Rscript

# Function that takes multiples gene count bed files and 
# merges them into on file (matrix). Intended for the mgs_augmented_assembly pipeline
# Author: Julien Tremblay - jtremblay514@gmail.com

options(stringsAsFactors = FALSE)
mergeTables <- function(infiles_htseq, infiles_cov, outfileRaw, outfileNormalized) {
  
  #infiles_htseq = "~/Projects/test_MT/395160.htseq,~/Projects/test_MT/395161.htseq"
  #infiles_cov = "~/Projects/test_MT/395160.genes.cov,~/Projects/test_MT/395161.genes.cov"
  
  files_htseq = unlist(strsplit(infiles_htseq, ","))
  #files_cov = unlist(strsplit(infiles_cov, ","))
  
  df = NULL
  dfNormalized = NULL
  
  # symbols will be the same in both .htseq and .cov
  tData = read.table(files_htseq[1], sep="\t", header=FALSE, comment.char="")
  symbol = as.character(tData[,1])
  
  df = cbind(df, symbol)
  df = data.frame(df)
  df = df[grep("__", df[,1], invert=TRUE),]
  df = data.frame(df)
  df = df[order(df),]
  df = data.frame(df)
  #dfNormalized = cbind(dfNormalized, symbol)
  dfNormalized = df
  
  colNames = c("cluster")
  
  for(i in 1:length(files_htseq)){
    print(paste0("Reading file ", files_htseq[i]))
    currFileHtseq = files_htseq[i];
    currFileCov = files_cov[i];
    currName = basename(currFileHtseq)
    currName = gsub(".htseq", "", currName)
    colNames = append(colNames, currName)
    
    # Use cov table to get gene length.
    currTableCov = read.table(currFileCov, sep="\t", header=FALSE, comment.char="")
    currTableCovSymbol = data.frame(currTableCov[,1])
    length = (currTableCov[,3] - currTableCov[,2]) / 1000 #To get gene lengths in kb.
    currDfCov = data.frame(currTableCov[,1])
    currDfCov = cbind(currDfCov, length)
    colnames(currDfCov) =  c("cluster", currName)
    
    # Then get total mapped htseq reads (dont forget to remove last rows starting with "__")
    currTableHtseq = read.table(currFileHtseq, sep="\t", header=FALSE, comment.char="")
    currTableHtseq2 = currTableHtseq[grep("__", currTableHtseq$V1, invert=TRUE), ]
    
    # then sort both tables to make sure gene_id order is the same in both tables.
    currTableHtseq3 = currTableHtseq2[order(currTableHtseq2[,1]),]
    currDfCov2 = currDfCov[order(currDfCov[,1]),]
    
    totalReads = sum(as.numeric(currTableHtseq3[,2]))
    rawCounts = currTableHtseq3[,2]
    
    #normalizedCount = (currTable[,4] / (length * 1000)) # report mapped reads per kb.
    #RPKM = (currTable[,4] / (length * 1000) / (totalReads / 1000000)) # report mapped reads per kb per million reads mapped (RPKM).
    RPKM = (1000000 * currTableHtseq3[,2]) / (currDfCov2[,2] * totalReads) # report mapped reads per kb per million reads mapped (RPKM).
    df = cbind(df, rawCounts)
    colnames(df)[ncol(df)] = currName
    dfNormalized = cbind(dfNormalized, RPKM)
    colnames(dfNormalized)[ncol(dfNormalized)] = currName
    #print(paste0("Total reads:", totalReads))
  }
  df = data.frame(df)
  dfNormalized = data.frame(dfNormalized)
  names(df) = colNames
  names(dfNormalized) = colNames
  df[is.na(df)] = 0
  dfNormalized[is.na(dfNormalized)] = 0
  
  write.table(df, file=outfileRaw, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  write.table(dfNormalized, file=outfileNormalized, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  print("Done mergingTables...")
  
} 

usage=function(errM) {
  cat("\nUsage : Rscript mergeHtseqAbundance.R [option] <Value>\n")
  cat("       -i        : List of htseq files separated by a comma\n")
  #cat("       -c        : List of bedtools cov files separated by a comma\n")
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
  } #else if (ARG[i] == "-c") {
    #infiles_cov=ARG[i+1]
  #}
    else if (ARG[i] == "-o") {
    outfile=ARG[i+1]
  } else if (ARG[i] == "-n") {
    outfileNormalized=ARG[i+1]
  }   
}

mergeTables(infiles_htseq, outfile, outfileNormalized)
