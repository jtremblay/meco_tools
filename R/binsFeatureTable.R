#!/usr/bin/env Rscript

# Function that generates a Feature table based on bins data.
binsFeatureTable <- function(summarized_bins_file, parsed_bins_file, abundance_file, outfile) {
   library(plyr)
   library(data.table)
   options(stringsAsFactors = FALSE)
   
   vColors = c('#00FF00','#FF8080','#FF00FF','#0000FF','#808282','#CCFFFF','#CCFFCC','#99CCFF','#CC99FF','#FFCC99','#3366FF','#33CCCC','#99CC00','#FF99CC','#FFCC00')
   vColors = c("#0000CD","#00FF00","#FF0000","#808080","#000000","#B22222","#40E0D0","#DAA520","#DDA0DD","#FF00FF","#00FFFF","#4682B4","#008000","#E6E6FA","#FF8C00","#80008B","#8FBC8F","#00BFFF","#FFFF00","#808000")
   
   tAbundance <- data.frame(fread(abundance_file, header=T, sep="\t"), check.names=FALSE)
   tParsedBins <- data.frame(fread(parsed_bins_file, header=F, sep="\t"), check.names=FALSE)
   tSummarizedBins <- data.frame(fread(summarized_bins_file, header=T, sep="\t"), check.names=FALSE)
   # The taxonomy is always the last 7 columns (Kingdom..Species), whatever the number of
   # quality columns before it (CheckM and CheckM2 report different columns).
   taxCols = (ncol(tSummarizedBins) - 6):ncol(tSummarizedBins)
   taxPrefixes = c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
 
   colnames(tAbundance)[1] = "V1" 
   print(head(tAbundance))
   #q();

   #print(head(tParsedBins))
   tSummarizedBins[[taxCols[1]]] = gsub("k__", "", tSummarizedBins[[taxCols[1]]])
   print("tSummarizedBins")
   print(head(tSummarizedBins))

   bins = unique(tParsedBins$V1)
   
   df = tAbundance[1,]
   df$V1 = NULL
   df = df[-1,]
  
   #print("head df")
   #print(head(df))

   for(bin in bins){
      # get contigs id
      #print("bin:")
      #print(bin)
      contigs = tParsedBins[tParsedBins$V1 == bin, 2]
      
      #print("contigs:")
      #print(contigs)

      # from abundance file. get rows corresponding to contig id and sum them up.
      contigs_values = tAbundance[tAbundance$V1 %in% contigs,]

      row.names(contigs_values) = contigs_values$V1
      contigs_values$V1 = NULL
      contigs_values2 = data.frame(t(colSums(contigs_values)), check.names=FALSE, row.names=bin)
      
      # Attach row to df
      df = rbind(df, contigs_values2)
   }
   #print("df after loop:")
   #print(df)
      
   # then add taxonomic lineages.
   for(i in seq_along(taxCols)){
      tSummarizedBins[[taxCols[i]]] = paste0(taxPrefixes[i], tSummarizedBins[[taxCols[i]]])
   }
  
   #print("tSummarizedBins")
   #print(tSummarizedBins)

   tLink = data.frame(tSummarizedBins[[1]])
   colnames(tLink) = c("#FEATURE_ID")
   tLink$Lineage = do.call(paste, c(tSummarizedBins[,taxCols], sep=";"))
   print("tLink")
   print(tLink)
   #tLink = unique(tLink)
   #head(tLink)
   
   # Merge df with tLink
   df2 = merge(df, tLink, by.x="row.names", by.y="#FEATURE_ID")
   colnames(df2)[1] = "#FEATURE_ID"
   colnames(df2)[ncol(df2)] = "taxonomy"
   
   write.table(df2, outfile, quote=FALSE, sep="\t", row.names=FALSE)
   
   print("[DEBUG] Completed generating Feature table...")
} 

usage=function(errM) {
   cat("\nUsage : Rscript binsFeatureTable.R [option] <Value>\n")
   cat("       -s        : summarized bins\n")
   cat("       -p        : parsed bins\n")
   cat("       -a        : contigs abundance\n")
   cat("       -o        : outfile (Feature table)\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 4) {
   usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
   if (ARG[i] == "-s") {
      summarized_bins=ARG[i+1]
   } else if (ARG[i] == "-p") {
      parsed_bins=ARG[i+1]
   }else if (ARG[i] == "-a") {
      abundance=ARG[i+1]
   }else if (ARG[i] == "-o") {
      outfile=mappingFile=ARG[i+1]
   }
}

binsFeatureTable(summarized_bins, parsed_bins, abundance, outfile)

