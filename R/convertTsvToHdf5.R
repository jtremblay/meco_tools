#!/usr/bin/env Rscript

# Convert gene abundance matrix in tsv to hdf5.
# Author: Julien Tremblay - jtremblay514@gmail.com

options(stringsAsFactors = FALSE)
convertTsvToHdf5 <- function(infile, outfile) {
  
    library(DelayedArray)
    library(HDF5Array)
    library(data.table)

    gene_abundance <- fread(infile)
    gene_names <- gene_abundance[[1]]
    counts_matrix <- as.matrix(gene_abundance[, -1])
    
    rownames(counts_matrix) <- gene_names

    writeHDF5Array(counts_matrix,
               filepath=outfile,
               name="counts",
               with.dimnames=TRUE,
               chunkdim=c(10000, 10),
               level=6)
 
    print("Done converting table...")
} 

usage=function(errM) {
  cat("\nUsage : Rscript generateCPMs.R [option] <Value>\n")
  cat("       -i        : abundance matrix of feature[i.e. genes, contigs]\n")
  cat("       -o        : outfile of in tsv\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 2) {
  usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-i") {
    infile = ARG[i+1]
  } else if (ARG[i] == "-o") {
    outfile = ARG[i+1]
  }
}
convertTsvToHdf5(infile, outfile)
