#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)

# Function to generate a OTU table based on blast results (against Greengenes).
# National Research Council Canada - Biomonitoring
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
generateOTUTableForPicrust <- function(infile_blast, infile_otu_table, outfile) {
   
   #infile_blast = "~/Projects/picrust/blast.out.besthit"
   #infile_otu_table = "~/Projects/picrust/otu_table.tsv"
   #outfile = "~/Projects/picrust/otu_table_picrust.tsv"
   
   tBlast = data.frame(fread(infile_blast, sep="\t", header=FALSE), check.names=FALSE)
   colnames(tBlast)[1] = "query"
   colnames(tBlast)[2] = "subject"
   otu_table = data.frame(fread(infile_otu_table, sep="\t", header=TRUE), check.names=FALSE)
   
   otu_table2 = merge(otu_table, tBlast[,c("query", "subject")], by.x="#OTU ID", by.y="query")
   otu_table2[["#OTU ID"]] = otu_table2$subject
   otu_table2$subject = NULL
   sum = rowSums(otu_table2[,2:(ncol(otu_table2)-1)])
   otu_table3 = otu_table2
   otu_table3$sum = sum
   otu_table3 = otu_table3[order(otu_table3$sum, decreasing=TRUE),]
   otu_table3$sum = NULL
   
   write.table(otu_table3, outfile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

}

usage=function(errM) {
  cat("\nUsage : Rscript mergeRawCounts.R [option] <Value>\n")
  cat("       -i        : List of files separated by a comma\n")
  cat("       -o        : outfile\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 2) {
  usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-b") {
    infile_blast=ARG[i+1]
  } else if (ARG[i] == "-t") {
    infile_otu_table=ARG[i+1]
  } else if (ARG[i] == "-o"){
    outfile=ARG[i+1]
  }
}

generateOTUTableForPicrust(infile_blast, infile_otu_table, outfile)

