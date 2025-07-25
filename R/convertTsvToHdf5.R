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

    dims <- dim(counts_matrix)
    n_rows <- dims[1]
    n_cols <- dims[2]

    # Choose a target chunk size (in number of elements). Adjust to your I/O needs.
    target_chunk_elems <- 100000  # ~400 KB for integer data

    # Compute chunk row size, default to full columns
    chunk_rows <- min(n_rows, max(1, floor(target_chunk_elems / n_cols)))

    # Construct chunkdim vector
    chunkdim_auto <- c(chunk_rows, n_cols)

    writeHDF5Array(counts_matrix,
               filepath=outfile,
               name="counts",
               with.dimnames=TRUE,
               chunkdim=chunkdim_auto,
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
