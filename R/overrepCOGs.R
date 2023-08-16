#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
library(data.table)
library(dplyr)

# National Research Council - Biomonitoring
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
overrepCOGs <- function(infile_abundance, infile_cogs, outfile) {

    cog = data.frame(fread(infile_cogs, header=FALSE), check.names=FALSE)
    abundance = data.frame(fread(infile_abundance, header=TRUE), check.names=FALSE)

    colnames(abundance)[1] = "cluster"

    print(head(abundance))
    print(head(cog))

    abundance = merge(abundance, cog[,c("V1", "V13")], by.x="cluster", by.y="V1")
    abundance$cluster = NULL


    abundance_COG = abundance %>%
        group_by(V13) %>% # means that the following summarize function will be done by the kegg_entry variable.
        dplyr::summarise_all(funs(sum)) %>%
        as.data.frame()

    colnames(abundance_COG)[1] = "COG"
    write.table(abundance_COG, outfile, sep="\t", quote=FALSE, row.names=FALSE)
             

}

usage=function(errM) {
  cat("\nUsage : Rscript overrepCOGs.R [option] <Value>\n")
  cat("       -i        : infile_abundance\n")
  cat("       -c        : infile_cogs\n")
  cat("       -o        : outfile\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 3) {
  usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-i") {
    infile_abundance=ARG[i+1]
  } else if (ARG[i] == "-c") {
    infile_cogs=ARG[i+1]
  }else if (ARG[i] == "-o") {
    outfile=ARG[i+1]
  }
}

overrepCOGs(infile_abundance, infile_cogs, outfile)
