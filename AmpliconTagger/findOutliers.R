#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# findOutliers.R.
# National Research Council Canada - Genomics and Microbiomes
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
findOutliers <- function(infile, prob, min_count) {
   options(stringsAsFactors = FALSE)
   library(outliers)
   library(data.table)
   #data_file = "~/Projects/Yergeau/FE/16S/countReport.tsv"
   #data_file = "~/Projects/Lallemand_AAD/countReport.tsv"
   #data_file = "~/Projects/PascalDrouin/16S/countReport.tsv"
   #infile = "~/build/misc_R_scripts/otu_table_filtered_bacteria.tsv"
   #prob=0.995
   #min_count=7000
   prob = as.numeric(prob)
   min_count = as.numeric(min_count)

   data = data.frame(fread(infile, sep="\t", header=TRUE, skip="#OTU ID"), check.names=FALSE)
   data = data.frame(count = colSums(data[,2:(ncol(data)-1)]))
   data$count = as.numeric(data$count)
   data$log = log(data$count + 1)
   data$name = row.names(data)
   data = data[,c("name", "count", "log")]
   row.names(data) = NULL
   data = data[order(data$log),]
   data$n = seq(1,nrow(data), 1)       
   
   cut = 0
   data2 = data[scores(data$log, type="z", prob=prob),]  
   print(data2)
   found_lower_tail = FALSE
   for(i in 1:nrow(data2)){
      if(i == data2[i,c("n")]){
         cut = i
         found_lower_tail = TRUE
      }else{
         print(paste0("i=: ", i, " is an index in the higher tail."))
      }
   }
   
   if(found_lower_tail == TRUE){
      data2 = data2[1:cut,] 
      threshold = data2[nrow(data2),]$count
      # Double check if lower tail contains enough reads to be kept - for instance above 10,000, just keep the lower value.
      for(i in 1:nrow(data2)){
         if(data2[i,]$count >= min_count){
            threshold = data2[i,]$count
            break;
         }
      }
   }else{
      threshold = data[1,c("count")]
   }
   
   print("The threshold that was found for rarefaction is :")
   print(paste0("THRESHOLD: ", threshold))
   print("The raw OTU table should be rarefied at this number of reads.")
   print("Here are the samples that will be filtered/left out if using this cutoff:")
   print(data2)
}

usage=function(errM) {
   cat("\nUsage : Rscript findOutliers.R [option] <value>\n")
   cat("       -i        : otu table\n")
   cat("       -p        : Z-score probability to determine which samples to exclude from OTU table\n")
   cat("       -n        : If the rarefaction value found by the Z-score is above the specified int value, the sample with the minimum number of reads above that specified int value will be picked as rarefaction value. --n <int> : If table is to be rarefied at <int> reads\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 3) {
   usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
   if (ARG[i] == "-i") {
      infile=ARG[i+1]
   } else if (ARG[i] == "-p") {
      prob=ARG[i+1]
   } else if (ARG[i] == "-m") {
      min_count=ARG[i+1]
   }
}

findOutliers(infile, prob, min_count) 



