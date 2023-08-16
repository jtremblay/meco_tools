#!/usr/bin/env Rscript

# Function that takes multiples feature tables and merges them into on file (matrix).
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca

library(data.table)
options(stringsAsFactors = FALSE)
mergeFeatureTablesIntoSingleFT <- function(infiles, outfile) {

  files = unlist(strsplit(infiles, ","))

  basename = gsub("_feature_table.tsv", "", files[1])
  tData = data.frame(fread(files[1], sep="\t", header=TRUE), check.names=FALSE)
  tData$`#FEATURE_ID` = paste0(basename, "_", tData$`#FEATURE_ID`)
  refTaxTable = data.frame(feature_id=tData$`#FEATURE_ID`, taxonomy=tData$taxonomy)
  tData$taxonomy = NULL

  for(i in 2:length(files)){
    print(paste0("Reading file ", files[i]))
    currFile = files[i];
    basename = gsub("_feature_table.tsv", "", files[i])
    currTable = data.frame(fread(currFile, sep="\t", header=TRUE), check.names=FALSE)
    currTable$`#FEATURE_ID` = paste0(basename, "_", currTable$`#FEATURE_ID`)
    refTaxTable = rbind(refTaxTable,
                        data.frame(feature_id=currTable$`#FEATURE_ID`, taxonomy=currTable$taxonomy)
    )
    currTable$taxonomy = NULL
    currMergedTable = merge(tData, currTable, by="#FEATURE_ID", all=T)
    tData = currMergedTable
  }

  df = NULL
  df = data.frame(tData, check.names=FALSE)
  df[is.na(df)] = 0
  df = merge(df, refTaxTable,by.x="#FEATURE_ID", by.y="feature_id", all=FALSE)

  write.table(df, file=outfileRaw, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

  print("Done merging feature tables...")
}

usage=function(errM) {
  cat("\nUsage : Rscript mergeFeatureTablesIntoSingleFT [option] <Value>\n")
  cat("       -i        : List of files separated by a comma\n")
  cat("       -o        : outfile raw counts\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 3) {
  usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-i") {
    infiles=ARG[i+1]
  } else if (ARG[i] == "-o") {
    outfile=ARG[i+1]
  }
}

mergeFeatureTablesIntoSingleFT(infiles, outfile)

