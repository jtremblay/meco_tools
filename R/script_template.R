#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# Tetra nucleotide frequency test.
# National Research Council - Biomonitoring
# Author: Julien Tremblay - jtremblay514@gmail.com
computeTetraNuclFreq <- function(infile, outfile) {

}

uiage=function(errM) {
  cat("\nUsage : Rscript mergeRawCounts.R [option] <Value>\n")
  cat("       -i        : List of files separated by a comma\n")
  cat("       -o        : outfile\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 2) {
  usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-i") {
    infiles=ARG[i+1]
  } else if (ARG[i] == "-o") {
    outfile=ARG[i+1]
  }
}

