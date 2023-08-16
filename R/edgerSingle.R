#!/usr/bin/env Rscript

# Function that takes a infile: (i.e. gene count table)  and normalize it using edgeR's RLE method.
# Julien Tremblay => julien.tremblay@nrc-cnrc.gc.ca
performEdgerSingle <- function(infile, outfile) {
   library(edgeR)
   library(tools)

   options(stringsAsFactors = FALSE)
   
   #infile = "~/Projects/metagenome_stats/merged_gene_abundance_50000.tsv"
   #outfile = "~/Projects/metagenome_stats/merged_gene_abundance_50000_normalized.tsv"
   
   print(paste0("[DEBUG] infile: ", infile))
   print(paste0("[DEBUG] outdir: ", outfile))
   
   #load data and work the gene count table to make it compatible with edgeR
   tData <- read.csv(file=infile, header=T, sep="\t", check.names=FALSE)
   colnames(tData)[1] = "Symbol"
   tData2 = tData[, 2:ncol(tData)]
   rownames(tData2) <- tData[ ,1]
   tData2$Symbol = NULL 
   
   x = tData2 + 1
   
   # EDGER
   y = edgeR::DGEList(counts = x)#, remove.zeros = TRUE)
   z = edgeR::calcNormFactors(y)
      
   x1 = x
   for(i in 1:ncol(x)){
      lib_size = z$samples[c(names(x)[i]),]$lib.size
      norm_fac = z$samples[c(names(x)[i]),]$norm.factors
      x1[,i] = x[,i] / (lib_size * norm_fac)
   }
      
   # DESEQ
   #cds = newCountDataSet(x, sampleConditions, featureData = taxADF)
   #cds = estimateSizeFactors(cds)
   #cds = estimateDispersions(cds, ...)
   
   # names
   #x1 = cbind(x1, tData[,ncol(tData)])
   #names(x1)[ncol(x1)] = "taxonomy"
   
   curr_row_names = rownames(x1)
   df = NULL
   df = cbind(df, curr_row_names)
   df = data.frame(df)
   df = cbind(df, x1[,1:ncol(x1)])
   names(df)[1] = "gene_id"
   
   #line="#OTU table normalized with edgeR RLE method"
   #write(line, file=outfile ,append=FALSE)
   write.table(df, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=FALSE)
   
   print("[DEBUG] Done Normalizing counts...")
} 

usage=function(errM) {
        cat("\nUsage : Rscript taxBarplotWithMappingFile.R [option] <Value>\n")
        cat("       -i        : infile\n")
        cat("       -o        : outfile\n")
        #cat("       -c        : cutoff\n")
        #cat("       -m        : design_file\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 2) {
	usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
	if (ARG[i] == "-i") {
		infile=ARG[i+1]
	} else if (ARG[i] == "-o") {
		outfile=ARG[i+1]
	#} else if (ARG[i] == "-c") {
	#  cutoff=ARG[i+1]
	}
}

performEdgerSingle(infile, outfile)
