#!/usr/bin/env Rscript

# generateCPMsH5. Normalize CPMs from HDF5.
# Author: Julien Tremblay - jtremblay514@gmail.com

options(stringsAsFactors = FALSE)
generateCPMsH5 <- function(infile, outfile){
  
    library(DelayedArray)
    library(HDF5Array)
    library(edgeR)
    library(data.table)

    # Assuming you have a large count matrix stored in a file (e.g., HDF5)
    # Replace with your actual file path and data structure
    message("Loading hdf5 array.")
    count_matrix <- HDF5Array(infile, "counts")

    message("creating DGEList.")
    # Create a DGEList object (from edgeR) with DelayedMatrix
    dge <- DGEList(counts=count_matrix)

    # Perform TMM normalization (will work with DelayedMatrix)
    message("Performing calcNormFactors.")
    dge <- calcNormFactors(dge, method="TMM")

    # generate CPMs
    message("Generating CPMs.")
    cpms_matrix_delayed <- cpm(dge, normalized.lib.sizes=TRUE, log=FALSE)
    message("Rounding.")
    #cpms_matrix_delayed <- round(cpms_matrix_delayed, digits=3)
    # Setup (before the loop)
    output_file <- "merged_gene_abundance_cpm.tsv"
    chunk_size <- 10000

    is_first_chunk <- TRUE

    for(i in seq(1, nrow(cpms_matrix_delayed), by=chunk_size)) {
        end_row <- min(i + chunk_size - 1, nrow(cpms_matrix_delayed))
    
        current_chunk <- as.matrix(cpms_matrix_delayed[i:end_row, ])
        current_chunk_rounded <- round(current_chunk, digits=3)
    
        # Convert to data.table for fwrite
        current_chunk_dt <- as.data.table(current_chunk_rounded, keep.rownames = TRUE)
        # The 'keep.rownames = TRUE' will convert row names to a column named "rn" by default.
        # You might want to rename this column to "GeneID" or similar later if needed.
    
        if(is_first_chunk){
            # Write the header and first batch of data
            fwrite(current_chunk_dt,
                   file=output_file,
                   append=FALSE, # Overwrite/create file for the first chunk
                   sep="\t",
                   row.names=FALSE, # Row names are now a column
                   col.names=TRUE, # Include column names (header)
                   quote=FALSE)
            is_first_chunk <- FALSE
        } else {
            # Append subsequent batches of data without the header
            fwrite(current_chunk_dt,
                   file=output_file,
                   append=TRUE,
                   sep="\t",
                   row.names=FALSE, # Row names are now a column
                   col.names=FALSE, # IMPORTANT: Do NOT include column names again
                   quote=FALSE)
        }
    
        rm(current_chunk, current_chunk_rounded, current_chunk_dt)
        gc()
    }

    # write to delayed array.
    #writeHDF5Array(counts_matrix,
    #           filepath=outfile_h5,
    #           name="counts",
    #           with.dimnames=TRUE)
 
    print("Done normalizing table...")
} 

usage=function(errM) {
  cat("\nUsage : Rscript generateCPMsH5.R [option] <Value>\n")
  cat("       -i        : abundance matrix in h5 format.\n")
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
generateCPMsH5(infile, outfile)
