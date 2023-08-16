#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
#library(dplyr)
library(data.table)

# National Research Council - Biomonitoring
# Author: Julien Tremblay - jtremblay514@gmail.com
# mapping file has to have a column UniqueTreatment
computeAverageFromMergedAbundance <- function(abundance_file, mapping_file, outfile) {
   # Read annotations and associate KO with each gene_id (if available)
   #root = "~/Projects/PascalDrouinMG/export/"
   #abundance_file = paste0(root, "/genes/merged_gene_abundance_cpm.tsv")
   #mapping_file = paste0(root, "/mapping_file.tsv")
   #outfile = paste0(root, "/genes/merged_gene_abundance_cpm_means.tsv")
   
   # Load abundance
   abundance = data.frame(fread(abundance_file, sep="\t", header=TRUE), check.names=FALSE)
   colnames(abundance)[1] = "gene_id"
   
   # load mapping file. 
   mapping = data.frame(fread(mapping_file, sep="\t", header=TRUE), check.names=FALSE)
   #mapping = mapping[mapping$UniqueTreatment %in% curr_treatments,]
   
   abundance = abundance[,colnames(abundance) %in% c("gene_id", mapping$`#SampleID`)]
   
   df = NULL
   treatments = unique(mapping$Treatment)
   for(i in 1:length(treatments)){
      treatment = treatments[i]
      curr_samples = as.character(mapping[mapping$Treatment == treatment,]$`#SampleID`)
      abundance2 = abundance[,c("gene_id", curr_samples)]
      if(ncol(abundance2) >= 3){
         abundance2$mean = round(rowMeans(abundance2[,2:ncol(abundance2)]), digits=2)
      }else{
         abundance2$mean = round(abundance2[,1:2,drop=FALSE], digits=2)
      }
      
      if(i == 1){
         df = abundance2[,c("gene_id", "mean"),drop=FALSE]
         colnames(df)[(i+1)] = treatment
      }else{
         df = cbind(df, abundance2$mean)
         colnames(df)[(i+1)] = treatment
      }
      #print(head(df))
      fwrite(df, outfile, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
      print("Done averaging treatments...")   
   }
}

uiage=function(errM) {
   cat("\nUsage : Rscript mergeRawCounts.R [option] <Value>\n")
   cat("       -i        : otu table\n")
   cat("       -n        : remove samples having less than n reads in total\n")
   cat("       -o        : outfile otu table\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 3) {
   usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
   if (ARG[i] == "-i") {
      infile      = ARG[i+1]

   } else if (ARG[i] == "-o") {
      outfile     = ARG[i+1]

   } else if (ARG[i] == "-m") {
      mapping_file= ARG[i+1]
   }
}
computeAverageFromMergedAbundance(infile, mapping_file, outfile)
