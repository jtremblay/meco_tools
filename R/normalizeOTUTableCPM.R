#!/usr/bin/env Rscript

# Function that takes a infile: (Qiime's OTU table in tsv)  and normalize it using edgeR's RLE method.
# Julien Tremblay => julien.tremblay@nrc-cnrc.gc.ca
performNormalizationOTUs <- function(infile, outfile, abundance_cpm) {
   library(edgeR)
   library(tools)
   library(data.table)

   options(stringsAsFactors = FALSE)
   
   #infile = "~/Projects/testRPlots/otu_table.tsv"
   #outfile = "~/Projects/testRPlots/otu_table_bacteriaArchaea_norm.tsv"
         
   print(paste0("[DEBUG] infile: ", infile))
   print(paste0("[DEBUG] outdir: ", outfile))
   
   #load data and work the gene count table to make it compatible with edgeR
   tData = fread(infile, skip="OTU ID", header=TRUE) # Means will start reading when find a row matching 'OTU ID'
   tData = data.frame(tData, check.names=FALSE)
   tData = tData[rowSums(tData[,2:(ncol(tData)-1)])!=0,]
   colnames(tData)[1] = "Symbol"
   #print(head(tData))
   #tData2 = tData[, 2:ncol(tData)-1]
   #rownames(tData2) <- tData[ ,1]
   #tData2$Symbol = NULL 
   
   tDataCPM = fread(abundance_cpm, header=TRUE)
   tDataCPM = data.frame(tDataCPM, check.names=FALSE)
   colnames(tDataCPM)[1] = "Symbol"
   rownames(tDataCPM) = tDataCPM$Symbol
   #tDataCPM$Symbol = NULL
   #print(head(tDataCPM))
   
   df = merge(tDataCPM, tData[,c("Symbol", "taxonomy")], by.x="Symbol", by.y="Symbol")
   colnames(df)[1] = "#OTU ID"
   #print(head(df))

   line="# normalized with edgeR (i.e. values = CPMs of complete abundance matrix(genes or contigs) )"
   write(line, file=outfile ,append=FALSE)
   write.table(df, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=TRUE)
   
   print("[DEBUG] Done Normalizing OTU table...")
} 

usage=function(errM) {
        cat("\nUsage : Rscript normalizeOTUTableCPM.R [option] <Value>\n")
        cat("       -i        : infile\n")
        cat("       -o        : outfile\n")
        cat("       -c        : abundance_cpm\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 3) {
	usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
	if (ARG[i] == "-i") {
		infile=ARG[i+1]
	} else if (ARG[i] == "-o") {
		outfile=ARG[i+1]
	} else if (ARG[i] == "-c") {
	  abundance_cpm=ARG[i+1]
	}
}

performNormalizationOTUs(infile, outfile, abundance_cpm)

