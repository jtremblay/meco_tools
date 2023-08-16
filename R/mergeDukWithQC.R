#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# Merge duk logs with qc mapping stats.
# National Research Council - Biomonitoring
# Author: Julien Tremblay - jtremblay514@gmail.com
mergeDukWithQC <- function(infile_duk, infile_qc, outfile) {

   library(data.table)
   
   duk = data.frame(fread(infile_duk, header=TRUE, showProgress=FALSE), check.names=FALSE)
   qc = data.frame(fread(infile_qc, header=TRUE, showProgress=FALSE), check.names=FALSE)

   tData = merge(qc, duk, by="sampleName")

   colnames(tData)[ncol(tData)] = "ContaminantRatio"
   colnames(tData)[ncol(tData)-1] = "ContaminantReads"

   new_order = c("sampleName", "rawFragments", "survivingFragments", "survivingFragments%", "survivingSingle", "totalReadsQCed",
              "ContaminantReads", "ContaminantRatio", "mapped", "mapped%", "properlyPaired", "properlyPaired%", "TOTAL")
   
   tData2 = tData[,new_order]


   df = NULL
   df = data.frame(mapped=as.numeric(gsub(",","",tData2$mapped)), 
                total=as.numeric(gsub(",","",tData$TOTAL)), 
                contaminants=as.numeric(gsub(",","",tData2$ContaminantReads)))

   tData$TOTAL = NULL
   df$total = df$total * 1


   tData2[["mapped%"]] = df$mapped_perc = ( df$mapped/(df$total - df$contaminants) ) * 100
   tData2[["mapped%"]] = paste0(round(tData2[["mapped%"]], 0),"%")

   tData2$ContaminantReads = format(tData2$ContaminantReads, big.mark=",", scientific=FALSE)
   write.table(tData2, file = outfile, sep="\t", quote=FALSE, row.names=FALSE)


}

usage=function(errM) {
  cat("\nUsage : Rscript mergeRawCounts.R [option] <Value>\n")
  cat("       -i        : contaminant duk summary file\n")
  cat("       -q        : qc reads / mapped reads file\n")
  cat("       -o        : outfile\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 3) {
  usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-i") {
    infile_duk=ARG[i+1]
  } else if (ARG[i] == "-q") {
    infile_qc=ARG[i+1]
  } else if (ARG[i] == "-o") {
    outfile=ARG[i+1]
  }
}

mergeDukWithQC(infile_duk, infile_qc, outfile)
