#!/usr/bin/env Rscript

# Function that takes a Feature table in input and writes that same Feature table sorted in abundance order (highest to lowest) in output.
# National Research Council Canada
# Genomics and Microbiomes
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
sortFeatureTable <- function(infile, outfile) {

    library(data.table)
    data = data.frame(fread(infile, skip="#FEATURE_ID", header=TRUE), check.names=FALSE)
    data = cbind(data, (rowSums(data[2:(ncol(data)-1)]))/(ncol(data)-2) )
    data = data[order(-data[, ncol(data)]),,drop=FALSE]
    data[,ncol(data)] = NULL
    write.table(data, outfile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
} 

usage=function(errM) {
        cat("\nUsage : Rscript sortFeatureTable.R [option] <Value>\n")
        cat("       -i        : Tree file\n")
        cat("       -o        : outfile\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 2) {
	usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
	if(ARG[i] == "-i") {
		infile=ARG[i+1]
	}else if (ARG[i] == "-o") {
		outfile=ARG[i+1]
	}
}

sortFeatureTable(infile, outfile)

