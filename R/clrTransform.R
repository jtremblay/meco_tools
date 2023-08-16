#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# Center Log-ration transform.
# National Research Council - Biomonitoring
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
clrTransform <- function(infile, outfile) {

    library(data.table)
    library(compositions)

    df = data.frame(fread(infile, header=TRUE, sep="\t"), check.names=FALSE)
    df$id = gsub("^(\\S+) .*", "\\1", df[,1])
    row.names(df) = df$id
    df$id = NULL
    df[,1] = NULL
    df_clr = as.data.frame(clr(df))
    write.table(df_clr, outfile, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)



}

usage=function(errM) {
  cat("\nUsage : Rscript clrTransform.R [option] <Value>\n")
  cat("       -i        : tetra nucleotide matrix to transform\n")
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

clrTransform(infile, outfile)
