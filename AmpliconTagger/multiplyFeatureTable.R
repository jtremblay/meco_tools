#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)
# Multiply Feature table.
# National Research Council Canada - Genomics and Microbiomes
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
multiplyFeatureTable <- function(infile, outfile, m) {
   
   #infile = "~/Projects/Lallemand_AAD/16S/export/feature_tables/feature_table_filtered.tsv"
   #outfile = "~/Projects/Lallemand_AAD/16S/export/feature_tables/feature_table_filtered_x1000.tsv"
   #m = 1000
   m = as.numeric(m)
   feature_table  = data.frame(fread(infile, sep="\t", header=TRUE, skip="#FEATURE_ID"), check.names=FALSE)
   feature_table2 = cbind(feature_table[,1], feature_table[,2:(ncol(feature_table)-1)] * m, feature_table[,ncol(feature_table)])
   colnames(feature_table2)[1] = "#FEATURE_ID"
   colnames(feature_table2)[ncol(feature_table2)] = "taxonomy"
   
   write.table(feature_table2, outfile, sep="\t", row.names=FALSE, quote=FALSE)
   
}

usage=function(errM) {
   cat("\nUsage : Rscript multiplyFeatureTable.R [option] <Value>\n")
   cat("       -i        : feature table\n")
   cat("       -m        : Multiplier\n")
   cat("       -o        : outfile feature table (multiplied)\n")
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
