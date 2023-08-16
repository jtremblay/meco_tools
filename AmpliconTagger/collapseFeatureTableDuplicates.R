#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

library(data.table)
library(dplyr)

# Aggregate duplicates of a Feature table (closed ref)
# National Research Council Canada
# Genomics and Microbiomes
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
collapseFeatureTableDuplicates <- function(infile, outfile) {
  
  #infile = "~/build/misc/feature_table_greengenes.tsv"
  #outfile = "~/build/misc/feature_table_greengenes_nodups.tsv"
  
  feature_table = data.frame(fread(infile, sep="\t", header=TRUE, skip="#FEATURE_ID"), check.names=FALSE)
  colnames(feature_table)[1] = "feature_id"
  
  link = feature_table[,c(1,(ncol(feature_table)))]
  link = unique(link)
  
  feature_table2 = feature_table %>%
    select(-one_of("taxonomy")) %>% 
      group_by(otu_id) %>% 
        summarize_all(sum) %>% 
          as.data.frame()
      
  feature_table2 = merge(feature_table2, link, by="otu_id")
  colnames(feature_table2)[1] = "#FEATURE_ID"
  feature_table2$total = rowSums(feature_table2[,2:(ncol(feature_table2)-1)])
  feature_table3 = feature_table2[order(-feature_table2$total),]
  feature_table3$total = NULL
  
  write.table(feature_table3, outfile, row.names=FALSE, sep="\t", quote=FALSE)
}

uiage=function(errM) {
  cat("\nUsage : Rscript collapseFeatureTableDuplicates.R [option] <Value>\n")
  cat("       -i        : Feature table file in tsv format\n")
  cat("       -o        : outfile\n")
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

collapseFeatureTableDuplicates(infile, outfile)
