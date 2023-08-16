#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(DropletUtils)
library(Matrix)

# Convert Single cell STAR output to hdf5 format.
# 
# Author: Julien Tremblay - jtremblay514@gmail.com
mtxToH5 <- function(infile_mtx, infile_barcodes, infile_features, gene_info_file, outfile, my_chemistry, my_genome) {
      gene_info = read.table(gene_info_file, sep="\t", header=F, skip=1)

      base = basename(infile_mtx)
      base = gsub("_matrix.mtx", "", base)

      # read files
      df = readMM(infile_mtx)
      barcodes = read.table(infile_barcodes, header=F)$V1
      features = read.table(infile_features, header=F)$V1
      # Here we assume that gene_info is the exact same order as features
      if(identical(features, gene_info$V1) == TRUE){
        features = gene_info$V2
      }else{
        stop("features and gene_info are not identical.")
      }

      row.names(df) = features
      colnames(df) = barcodes

      # Here let's write the files in .h5 format for better portability.
      write10xCounts(
        outfile,
        df,
        #barcodes = colnames(x),
        #gene.id = rownames(x),
        #gene.symbol = gene.id,
        gene.type = "Gene Expression",
        overwrite = TRUE,
        type = "HDF5",
        genome = my_genome,
        version = c("3"),
        chemistry = my_chemistry,
        original.gem.groups = 1L,
        library.ids = "custom"
      )
}

usage=function(errM) {
  cat("\nUsage : Rscript mergeRawCounts.R [option] <Value>\n")
  cat("       -m        : .mtx <file>\n")
  cat("       -b        : barcodes <file>\n")
  cat("       -f        : features <file>\n")
  cat("       -c        : chemistry info <string>\n")
  cat("       -g        : genome info <string>\n")
  cat("       -o        : .h5 outfile <file>\n")
  cat("       -k        : gene info <file>\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 7) {
  usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-m") {
    mtx=ARG[i+1]
  } else if (ARG[i] == "-b") {
    barcodes=ARG[i+1]
  } else if (ARG[i] == "-f") {
    features=ARG[i+1]
  } else if (ARG[i] == "-c") {
    chemistry=ARG[i+1]
  } else if (ARG[i] == "-g") {
    genome=ARG[i+1]
  } else if (ARG[i] == "-o") {
    outfile=ARG[i+1]
  } else if (ARG[i] == "-k") {
    gene_info=ARG[i+1]
  }
}
mtxToH5(mtx, barcodes, features, gene_info, outfile, chemistry, genome)
