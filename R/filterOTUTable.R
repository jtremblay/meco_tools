#!/usr/bin/env Rscript

# Function that takes a infile: (Qiime's tax summary spreadsheets) and an outdir where to
# To write filtered otu table by abundance.
# Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
generatePlots <- function(infile, outdir, abundance_cutoff) {
   library(ggplot2)
   library(scales)
   library(reshape2)
   library(grid)
   library(data.table)
   #library(igraph)
   library(plyr)
   library(tools)
   options(stringsAsFactors = FALSE)
   
   
   #infile = "~/Projects/networks/otu_table_filtered_bacteriaArchaea_normalized.tsv"
   #outdir = "/home/jtrembla/Projects/networks/"
   #prefix = basename(file_path_sans_ext(infile))
   #abundance_cutoff = 1 #between 1 and 3 idealy...
  
   print(paste0("[DEBUG] infile: ",infile))
   print(paste0("[DEBUG] outdir: ",outdir))
   print(paste0("[DEBUG] prefix: ",prefix))
   #print(paste0("[DEBUG] prefix: ",mappingFile))
   
   #outfile1 <- paste0(outdir, prefix, "_", ".pdf")
   #outfile2 <- paste0(outdir, prefix, "_", ".jpeg")

   # For col annotations
   vColors = c(
     "#0000CD", "#00FF00", "#FF0000", "#808080", "#000000", "#B22222", "#40E0D0", "#DAA520", "#DDA0DD", "#FF00FF",
     "#00FFFF", "#4682B4", "#008000", "#E6E6FA", "#FF8C00", "#80008B", "#8FBC8F", "#00BFFF", "#FFFF00", "#808000"
   )
   
   # For row annotations
   vColors2 = c(
      "#FFCCCC", "#FFE5CC", "#FFFFCC", "#E5FFCC", "#CCFFCC", "#CCFFE5", "#CCFFFF", "#CCE5FF", "#CCCCFF", "#E5CCFF", "#FFCCFF", "#FFCCE5", "#FFFFFF",
      "#990000", "#666600", "#006666", "#330066", "#A0A0A0", "#99004C"
   ) 
   
   # Read OTU table (TSV format)
   tData = fread(infile, skip="OTU ID") # Means will start reading when find a row matching 'OTU ID'
   tData = data.frame(tData, check.names=FALSE)
   colnames(tData)[1] = "Symbol"
   tData$Symbol = as.character(tData$Symbol)
   tLineages = data.frame(tData$Symbol, tData[,ncol(tData)])
   colnames(tLineages) = c("Symbol", "taxonomy")
   tData2 = tData[, 1:ncol(tData)-1]
   rownames(tData2) <- tData$Symbol
   tData2$Symbol = NULL
   
   # Filter using various cutoff.
   cutoffs = c(0.0001, 
               0.00025, 
               0.0005, 
               0.00075
   ) # Will keep OTUs being abundant more than these percentages
   tData2_perc = data.frame(prop.table(as.matrix(tData2), 2))
      
   for(cutoff in cutoffs){
      # filter by cutoff
      tData2$filt<-apply(tData2_perc, 1, function(x) sum(x > cutoff)) #you can change the greater than to less than if you want to invert the count.
      df<-tData2[tData2$filt>=abundance_cutoff,] #the 2 is madxe up by me for the case of wanting 2 or more columns that are .005 or greater.  Change the 2 for your needs
      #tData2$filt<-NULL #deleting dummy columns
      df$filt<-NULL
      df$Symbol = row.names(df)
      # Put back taxonomy lineages.
      df2 = merge(df, tLineages, by="Symbol")
      names(df2)[1] = "#OTU ID"
      df2 = arrange(df2, -df2[,2])
   
      outfile = paste0(outdir,"/", prefix, "_filtered_", abundance_cutoff, "_X_", cutoff, ".tsv")
         
      line = paste0("# filtered by ", abundance_cutoff, "_X_", cutoff, " abundance cutoff")
      write(line, file=outfile ,append=FALSE)
      write.table(df2, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=TRUE)
   }
	
	print("[DEBUG] Printed plots...")
} 

usage=function(errM) {
        cat("\nUsage : Rscript taxBarplotWithMappingFile.R [option] <Value>\n")
        cat("       -i        : infile\n")
        cat("       -o        : outdir\n")
        cat("       -c        : cutoff\n")
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
		outdir=ARG[i+1]
	}else if (ARG[i] == "-p") {
		prefix=ARG[i+1]
	}else if (ARG[i] == "-c") {
		mappingFile=ARG[i+1]
	}
}

generatePlots(infile, outdir, prefix, mappingFile)

