#!/usr/bin/env Rscript

# Function that takes multiples taxonomy tables and merges them.
# merges them into on file (matrix). Intended for the mgs_augmented_assembly pipeline
# Author: Julien Tremblay - jtremblay514@gmail.com

library(data.table)
options(stringsAsFactors = FALSE)
mergeFeatureTablesIntoTaxTables <- function(infiles, outfile) {

  files = unlist(strsplit(infiles, ","))

  tData = data.frame(fread(files[1], sep="\t", header=TRUE), check.names=FALSE)
  #tData$`#FEATURE_ID` = NULL
  tData[,1] = NULL

  for(i in 2:length(files)){
    taxonomy = NULL
    print(paste0("Reading file ", files[i]))
    currFile = files[i];
    currTable = data.frame(fread(currFile, sep="\t", header=TRUE), check.names=FALSE)
    #currTable$`#FEATURE_ID` = NULL
    currTable[,1] = NULL
    currMergedTable = merge(tData, currTable, by="taxonomy", all=T)
    tData = currMergedTable
  }

  df = NULL
  df = data.frame(tData, check.names=FALSE)
  #names(df) = colNames
  df[is.na(df)] = 0
  # Then, the magic happens here with dplyr.
  df2 = df %>%
    group_by(taxonomy) %>% # means that the following summarize function will be done by the kegg_entry variable.
    summarise_all(list(sum)) %>%
    as.data.frame()

  colnames(df2)[1] = "Taxon"

  write.table(df2, file=outfileRaw, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

  print("Done merging taxonomy tables...")

}

usage=function(errM) {
  cat("\nUsage : Rscript mergeFeatureTablesIntoTaxTables.R [option] <Value>\n")
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
mergeFeatureTablesIntoTaxTables(infiles, outfile)
