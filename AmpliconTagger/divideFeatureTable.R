#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)
# Multiply Feature table.
# National Research Council Canada - Genomics and Microbiomes
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
multiplyFeatureTable <- function(infile, outfile, m) {
   
   #infile = "~/Projects/Lallemand_AAD/16S/export/otu_tables/otu_table_filtered.tsv"
   #outfile = "~/Projects/Lallemand_AAD/16S/export/otu_tables/otu_table_filtered_x1000.tsv"
   #m = 1000
   m = as.numeric(m)
   otu_table  = data.frame(fread(infile, sep="\t", header=TRUE, skip="#FEATURE_ID"), check.names=FALSE)
   otu_table2 = cbind(otu_table[,1], otu_table[,2:(ncol(otu_table)-1)] / m, otu_table[,ncol(otu_table)])
   colnames(otu_table2)[1] = "#FEATURE_ID"
   colnames(otu_table2)[ncol(otu_table2)] = "taxonomy"
   
   write.table(otu_table2, outfile, sep="\t", row.names=FALSE, quote=FALSE)
   
}

usage=function(errM) {
   cat("\nUsage : Rscript multiplyFeatureTable.R [option] <Value>\n")
   cat("       -i        : otu table\n")
   cat("       -m        : Multiplier\n")
   cat("       -o        : outfile otu table (multiplied)\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 3) {
   usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
   if (ARG[i] == "-i") {
      infile=ARG[i+1]
   } else if (ARG[i] == "-o") {
      outfile=ARG[i+1]
   } else if (ARG[i] == "-m") {
      m=ARG[i+1]
   }
}

multiplyFeatureTable(infile, outfile, m)
