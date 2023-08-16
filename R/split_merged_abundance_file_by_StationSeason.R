#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# splitGeneAbundanceMatrix.
# National Research Council - Biomonitoring
# Author: Julien Tremblay - jtremblay514@gmail.com
splitGeneAbundanceMatrix <- function(abundance_file, mapping_file, outdir) {
   options(stringsAsFactors = FALSE) 
   library(data.table)
   
   #abundance_file = "~/Projects/DFO/microcosms/export/genes/merged_gene_abundance_10000.tsv"
   #mapping_file = "~/Projects/DFO/microcosms/mapping_file.tsv"
   #outdir = "~/Projects/DFO/microcosms/export/genes"
   
   abundance = data.frame(fread(abundance_file, sep="\t", header=TRUE), check.names=FALSE)
   mapping = data.frame(fread(mapping_file, sep="\t", header=TRUE), check.names=FALSE)
   
   mapping$StationSeason = paste0(mapping$Station, "", mapping$Season)
   stationseasons = unique(mapping$StationSeason)
   
   for(stationseason in stationseasons){
      curr_samples = mapping[mapping$StationSeason == stationseason,]$`#SampleID`
      curr_abundance = abundance[,c("cluster", curr_samples)]
      fwrite(curr_abundance, file=paste0(outdir, "/merged_gene_abundance_", stationseason, ".tsv"), sep="\t")
   }
   
}

usage=function(errM) {
   cat("\nUsage : Rscript splitGeneAbundanceMatrix [option] <Value>\n")
   cat("       -i        : Gene abundance matrix - first column is labeled 'cluster' \n")
   cat("       -o        : outdir\n")
   cat("       -m        : mapping_file\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 3) {
   usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
   if (ARG[i] == "-i") {
      infile=ARG[i+1]
   }else if(ARG[i] == "-m"){ 
      mapping_file=ARG[i+1]
   }else if (ARG[i] == "-o") {
      outdir=ARG[i+1]
   }
}
splitGeneAbundanceMatrix(infile, mapping_file, outdir)