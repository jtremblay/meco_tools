#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)

# National Research Council - Biomonitoring
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
roundOTUTable <- function(infile, outfile) {

    data = data.frame(fread(infile, skip="#OTU ID"), check.names=FALSE)
    data[,2:(ncol(data)-1)] = round(data[,2:(ncol(data)-1)], digits=0)
    write.table(data, outfile, row.names=FALSE, quote=FALSE, sep="\t")

}

usage=function(errM) {
  cat("\nUsage : Rscript roundOTUTable.R [option] <Value>\n")
  cat("       -i        : infile\n")
  cat("       -o        : outfile\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 2) {
  usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-i") {
    infile=ARG[i+1]
  } else if (ARG[i] == "-o") {
    outfile=ARG[i+1]
  }
}

roundOTUTable(infile, outfile)
