#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# Tetra nucleotide frequency test.
# National Research Council - Biomonitoring
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
splitAbundanceForMaxbin <- function(infile, outdir) {
    library(data.table)

    #df = data.frame(fread(infile), check.names=FALSE)
    df = fread(infile)
    df$contigLen = NULL
    df$totalAvgDepth = NULL

    dropcols = grep(".bam-var$", colnames(df))
    df[, (dropcols) := NULL]   
    df = data.frame(df, check.names=FALSE)
    row.names(df) = df$contigName
    df$contigName = NULL

    print(head(df))
    for(i in 1:ncol(df)){
        df2 = df[,i, drop=FALSE]
        curr_colname = colnames(df2)
        curr_colname = gsub(".bam", "", curr_colname)
        print("df2:")
        print(head(df2))
        write.table(df2, paste0(outdir, "/", curr_colname, ".txt"), row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t")
    }
}

uiage=function(errM) {
  cat("\nUsage : Rscript splitAbundanceForMaxbin.R [option] <Value>\n")
  cat("       -i        : List of files separated by a comma\n")
  cat("       -o        : outdir\n")
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
    outdir=ARG[i+1]
  }
}

splitAbundanceForMaxbin(infile, outdir)
