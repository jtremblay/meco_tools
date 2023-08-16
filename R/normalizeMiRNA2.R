#!/usr/bin/env Rscript

# Takes a miRNA raw abundance matrix as input and generates
# rpkm and filtered tables.
# Author: Julien Tremblay - jtremblay514@gmail.com

options(stringsAsFactors = FALSE)
normalizeMiRNA2 <- function(infile, length_file, outfileRpkm) {

  library(edgeR)
  library(data.table)

  df = data.frame(fread(infile, sep="\t", header=TRUE, showProgress=FALSE), check.names=FALSE)
  df_length = data.frame(fread(length_file, sep="\t", header=FALSE, showProgress=FALSE), check.names=FALSE)

  print(head(df))
  print(head(df_length))

  df = merge(df, df_length, by.x="#CLUSTER", by.y="V1")
  length_vec = df$V2

  row.names(df) = df[,1]
  df[,1] = NULL
  df$V2 = NULL
  annotations_df = df[,c("mirbase_hairpin", "mirbase_mature")]
  annotations_df$`#CLUSTER` = row.names(annotations_df)
  df$mirbase_hairpin = NULL
  df$mirbase_mature = NULL
  colnames(df) = gsub("_mapped.fastq.gz", "", colnames(df))

  y <- DGEList(df, remove.zeros=FALSE)
  y <- calcNormFactors(y, method="TMM")
  rpkms = rpkm(y, gene.length = length_vec)
  rpkms = round(rpkms, digits=3)
  rpkms = data.frame(rpkms)

  orig_colnames = colnames(rpkms)
  rpkms$`#CLUSTER` = row.names(rpkms)
  rpkms = rpkms[,c("#CLUSTER", orig_colnames)]

  rpkms = merge(rpkms, annotations_df, by.x="#CLUSTER", by.y="#CLUSTER")

  fwrite(as.data.frame(rpkms), outfileRpkm, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

  
  print("Done normalizing tables...")
} 

usage=function(errM) {
  cat("\nUsage : Rscript normalizeMiRNA2.R [option] <Value>\n")
  cat("       -i        : infile abundance matrix\n")
  cat("       -l        : length file (for computing rpkms\n")
  cat("       -o        : outfile rpkm counts\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 5) {
  usage("missing arguments")
}

skipNorm = NULL
## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-i") {
    infile = ARG[i+1]
  } else if (ARG[i] == "-o") {
    outfileRpkm = ARG[i+1]
  } else if (ARG[i] == "-l") {
    length_file = ARG[i+1]
  }
}

normalizeMiRNA2(infile, length_file, outfileRpkm)
