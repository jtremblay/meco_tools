#!/usr/bin/env Rscript

# Takes a miRNA raw abundance matrix as input and generates
# rpkm and filtered tables.
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca

options(stringsAsFactors = FALSE)
mergeMiRNATables <- function(infile_abundance, infile_mature, outfile) {

  library(data.table)


  df = data.frame(fread(infile_abundance, sep="\t", header=TRUE, showProgress=FALSE), check.names=FALSE)
  #hairpin = data.frame(fread(infile_hairpin, sep="\t", header=FALSE, showProgress=FALSE), check.names=FALSE)
  mature = data.frame(fread(infile_mature, sep="\t", header=FALSE, showProgress=FALSE), check.names=FALSE)

  gene_ids = unique(c(mature$V1))
  df = df[df$`#CLUSTER` %in% gene_ids,]

  #df = merge(df, hairpin[,c("V1", "V2")], by.x="#CLUSTER", by.y="V1", all=TRUE)
  #colnames(df)[ncol(df)] = "mirbase_hairpin"
  df = merge(df, mature[,c("V1", "V2")], by.x="#CLUSTER", by.y="V1", all=TRUE)
  colnames(df)[ncol(df)] = "mirbase_mature"
  print(head(df))

  fwrite(as.data.frame(df), outfile, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

  # Then filter based on cutoff
  
  print("Done merging and normalizing tables...")
} 

usage=function(errM) {
  cat("\nUsage : Rscript mergeMiRNATables.R [option] <Value>\n")
  cat("       -a        : infile abundance matrix\n")
  #cat("       -h        : hairpin mirna\n")
  cat("       -m        : mature mirna\n")
  cat("       -o        : outfile\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 3) {
  usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-a") {
    infile_abundance = ARG[i+1]
  } else if (ARG[i] == "-o") {
    outfile = ARG[i+1]
  #} else if (ARG[i] == "-h") {
    #infile_hairpin = ARG[i+1]
  } else if (ARG[i] == "-m") {
    infile_mature = ARG[i+1]
  }
}

mergeMiRNATables(infile_abundance, infile_mature, outfile)
