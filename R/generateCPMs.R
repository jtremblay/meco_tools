#!/usr/bin/env Rscript

# Function that takes a reads-based count matrix and writes in output
# a Count-Per-Million matrix.
# Author: Julien Tremblay - jtremblay514@gmail.com

options(stringsAsFactors = FALSE)
generateCPMs <- function(infile, outfile, skip_norm_factors) {
  
  #infile = "/project/6049199/projects/mock_community_shotgunMG/contig_abundance/merged_contig_abundance.tsv"
  #outfile = "/project/6049199/projects/mock_community_shotgunMG/contig_abundance/merged_contig_abundance_cpm.tsv"
  #skip_norm_factors = "FALSE"
    
  library(edgeR)
  library(data.table)

  skip_norm_factors = as.logical(skip_norm_factors)
  message("skip_norm_factors: ", skip_norm_factors)

  df = data.frame(fread(infile, header=TRUE, sep="\t"), check.names=FALSE)

  row.names(df) = df[,1]
  df$cluster = NULL
  df$contig_id = NULL
  df$gene_id = NULL
  df$feature_id = NULL
  
  for(i in 1:ncol(df)){
     df[,i] = as.numeric(df[,i])
  }

  ## edgeR. Load table, remove low cpm. 
  df = df[which(rowSums(df) > 0),]
  print("Running normalization...")
  if(isFALSE(skip_norm_factors)){
      y = DGEList(df, remove.zeros=TRUE)
      y = calcNormFactors(y, method="TMM")
      cpms = cpm(y)
      cpms = round(cpms, digits=3)
  }else{
      cpms = data.frame(apply((df+0), 2, function(x) (round(((x/(sum(x)+0))*1000000), digits=3))))
      # gives exactly the same value as edgeR's cpm() function
      # But we'll skip this for now.
      #for (i in seq_along(cpms)){
      #  set(cpms, j=i, value=cpms[[i]] / norm_factors$norm.factors[i])
      #}
  }
  fwrite(as.data.frame(cpms), outfile, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
  
  print("Done normalizing table...")
} 

usage=function(errM) {
  cat("\nUsage : Rscript generateCPMs.R [option] <Value>\n")
  cat("       -i        : abundance matrix of feature[i.e. genes, contigs]\n")
  cat("       -o        : outfile of CPMs\n")
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
  } else if (ARG[i] == "-n") {
    skip_norm_factors = ARG[i+1]
  }
}

generateCPMs(infile, outfile, skip_norm_factors)
