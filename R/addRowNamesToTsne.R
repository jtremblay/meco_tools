#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# Center Log-ration transform.
# National Research Council - Biomonitoring
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
addRowNamesToTsne <- function(infile_tsne, infile_clr, outfile) {

    library(data.table)

    tsne = data.frame(fread(infile_tsne, header=FALSE, sep="\t"), check.names=FALSE)
    clr = data.frame(fread(infile_clr, header=TRUE, sep="\t", select=c(1,2)), check.names=FALSE)
    
    print("head(tsne)")
    print(head(tsne))

    print("head(clr)")
    print(head(clr))
   
    row.names(tsne) = clr$V1
    colnames(tsne) = c("D1", "D2")

    write.table(tsne, outfile, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
}

usage=function(errM) {
  cat("\nUsage : Rscript addRowNamesToTsne.R [option] <Value>\n")
  cat("       -i        : tetra nucleotide matrix to transform\n")
  cat("       -o        : outfile\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 2) {
  usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-t") {
    infile_tsne=ARG[i+1]
  } else if (ARG[i] == "-c") {
    infile_clr=ARG[i+1]
  } else if(ARG[i] == "-o") {
    outfile=ARG[i+1]
  }
}

addRowNamesToTsne(infile_tsne, infile_clr, outfile)
