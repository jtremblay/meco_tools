#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)
# Divide Feature table.
# Author: Julien Tremblay - jtremblay514@gmail.com
percentageFeatureTable <- function(infile, outfile, m) {
   
    #infile = "~/Projects/Lallemand_AAD/16S/export/otu_tables/otu_table_filtered.tsv"
    #outfile = "~/Projects/Lallemand_AAD/16S/export/otu_tables/otu_table_filtered_x1000.tsv"
    # Here, we assume we have a header.
    df = data.frame(fread(infile, sep="\t", header=T), check.names=F)
    head(df)
    
    metadata = df[,c(1,ncol(df))]

    df[,c(1,ncol(df))] = NULL
    df = sweep(df, 2, colSums(df), "/")
    df = cbind(metadata[,1], df, metadata[,2])
    colnames(df)[1] = "#FEATURE_ID"
    colnames(df)[ncol(df)] = "taxonomy"
    write.table(df, outfile, sep="\t", quote=F, row.names=F)
}

usage=function(errM) {
   cat("\nUsage : Rscript multiplyFeatureTable.R [option] <Value>\n")
   cat("       -i        : otu table\n")
   cat("       -o        : outfile otu table (each cell divided by the sum of their respective column)\n")
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

percentageFeatureTable(infile, outfile)
