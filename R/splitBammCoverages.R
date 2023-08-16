#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# Split bamm parse output in individual files for maxbin compatibility.
# National Research Council - Biomonitoring
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
splitBammCoverage <- function(infile, outfile) {
  library(data.table)
  library(tools)

  #infile = "~/Projects/R/tmp.txt"
  #outdir = "~/Projects/R"
  
  tData = data.frame(fread(infile), check.names=FALSE)
  
  for(i in 3:ncol(tData)){
    prefix = basename(file_path_sans_ext(colnames(tData)[i]))
    filename = paste0(outdir, "/", prefix, ".txt")
    tData2 = data.frame(tData[,1], tData[,i])
    names(tData2) = c("contig", prefix)
    write.table(tData2, file=filename, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
    print(filename)
  }
  
}

uiage=function(errM) {
  cat("\nUsage : Rscript splitBammCoverage.R [option] <Value>\n")
  cat("       -i        : output of bamm parse\n")
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

splitBammCoverage(infile, outdir)
