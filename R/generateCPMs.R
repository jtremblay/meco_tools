#!/usr/bin/env Rscript

# Function that takes a reads-based count matrix and writes in output
# a Count-Per-Million matrix.
# Author: Julien Tremblay - jtremblay514@gmail.com

options(stringsAsFactors = FALSE)
generateCPMs <- function(infile, outfile) {

  library(edgeR)
  library(data.table)

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
  #tData = tData[which(rowSums(tData) > 0),]
  print("Running normalization...")
  y <- DGEList(df, remove.zeros=FALSE)
  y <- calcNormFactors(y, method="TMM")
  cpms = cpm(y)
  cpms = round(cpms, digits=3)
  #keep <- rowSums(cpm(y)>1) >= 2
  #y <- y[keep, , keep.lib.sizes=FALSE]
  #cpms = apply(df,2, function(x) (round( ((x/(sum(x)+1))*1000000), digits=3)))
  #cpms = round(cpm(y1), digits=3)
   
  #cpms = round(cpms, digits=3)
  #write.table(cpms, outfileCpm, quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)
  fwrite(as.data.frame(cpms), outfile, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
  
  print("Done merging and normalizing tables...")
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
  }
}

generateCPMs(infile, outfile)
