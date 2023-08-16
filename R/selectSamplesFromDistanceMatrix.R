#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# extract samples from a distance matrix.
# National Research Council Canada - Genomics and Microbiomes
# Author: Julien Tremblay - jtremblay514@gmail.com
library(data.table)
parseDistanceMatrix <- function(infile, outfile, samples) {
   
   #infile =  "~/build/misc_R_scripts/weighted_unifrac_otu_table_final_rarefied.txt"
   #outfile = "~/build/misc_R_scripts/weighted_unifrac_otu_table_final_rarefied_parsed.txt"
   #samples = "~/build/misc_R_scripts/samples_home_microbiome.txt"
   
   df = data.frame(fread(infile), check.names=FALSE)
   row.names(df) = df$V1
   df$V1 = NULL
   sample_names = data.frame(fread(samples, header=FALSE), check.names=FALSE)$V1
   
   df = df[row.names(df) %in% sample_names,]
   df = df[,colnames(df) %in% sample_names]
   
   # Then just to be sure, order them in the same order.
   missing_samples = sample_names[!sample_names %in% row.names(df)]
   sample_names = sample_names[sample_names %in% row.names(df)]
   
   df = df[,sample_names]
   df = df[sample_names,]
   
   curr_row_names = row.names(df)
   df$V1 = row.names(df)
   df = df[,c("V1", curr_row_names)]
   colnames(df)[1] = ""
   
   write.table(df, outfile, sep="\t", row.names=FALSE, quote=FALSE)
   
}

usage=function(errM) {
   cat("\nUsage : Rscript parseDistanceMatrix.R [option] <Value>\n")
   cat("       -i        : distance matrix\n")
   cat("       -s        : list of samples names in a text file. One sample name per row\n")
   cat("       -o        : outfile distance matrix\n")
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
   } else if (ARG[i] == "-s") {
      samples=ARG[i+1]
   }
}

parseDistanceMatrix(infile, outfile, samples)
