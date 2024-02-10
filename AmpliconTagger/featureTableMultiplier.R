#!/usr/bin/env Rscript

library(data.table)

options(stringsAsFactors = FALSE)

# Convert Single cell STAR output to hdf5 format.
# 
# Author: Julien Tremblay - jtremblay514@gmail.com
featureTableMultiplier <- function(infile, outfile, multiplier) {
    multiplier = as.numeric(multiplier)
    df = data.frame(fread(infile, sep="\t"), check.names=F)

    for(i in 2:(ncol(df)-1)){
        df[,i] = df[,i] * multiplier
    }
    write.table(df, outfile, sep="\t", quote=F, row.names=F)
}

usage=function(errM) {
  cat("\nUsage : Rscript featureTableX1000.R [option] <Value>\n")
  cat("       -i        : feature table <file>\n")
  cat("       -o        : feature table multiplied by multiplier  <file>\n")
  cat("       -n        : multiplier. default = 1000  <file>\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 3) {
  usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-i") {
    infile=ARG[i+1]
  } else if (ARG[i] == "-o") {
    outfile=ARG[i+1]
  } else if (ARG[i] == "-n") {
    multiplier=ARG[i+1]
  } 
}
featureTableMultiplier(infile, outfile, multiplier)
