#!/usr/bin/env Rscript

# Function that takes a infile: (Qiime's Feature table in tsv)  and normalize it using edgeR's RLE method.
# National Research Council Canada
# Genomics and Microbiomes
# Julien Tremblay => julien.tremblay@nrc-cnrc.gc.ca
performNormalizationFeatures <- function(infile, outfile, cutoff) {
   library(edgeR)
   library(tools)
   library(data.table)

   options(stringsAsFactors = FALSE)
   
   #infile = "~/Projects/test_edger_Features/feature_table_filtered_bacteriaArchaea.tsv"
   #outfile = "~/Projects/test_edger_Features/feature_table_filtered_bacteriaArchaea_norm.tsv"
   #outfile = "~/Projects/testRPlots/feature_table_bacteriaArchaea_norm.tsv"
   #cutoff = as.numeric(250)
   #multiplier = as.numeric(100000000)
   
   cutoff = as.numeric(cutoff)
   #multiplier = as.numeric(multiplier)
         
   print(paste0("[DEBUG] infile: ", infile))
   print(paste0("[DEBUG] outfile: ", outfile))
   print(paste0("[DEBUG] cutoff: ", cutoff))
   #print(paste0("[DEBUG] multiplier: ", multiplier))

   #load data and work the gene count table to make it compatible with edgeR
   #tData <- read.csv(file=infile, header=T, skip=1, sep="\t", check.names=FALSE)
   tData = data.frame(fread(infile, header=TRUE, skip="#FEATURE_ID"), check.names=FALSE) # Means will start reading when find a row matching 'FEATURE_ID'
   #tData = data.frame(tData, check.names=FALSE)
   print(head(tData))
   tData = tData[rowSums(tData[,2:(ncol(tData)-1)])!=0,]
   colnames(tData)[1] = "Symbol"
   ## Taxonomy equivalency
   tTaxonomy = tData[,c("Symbol", "taxonomy")]
   ##
   tData2 = tData[, 2:ncol(tData)-1]
   rownames(tData2) <- tData[ ,1]
   tData2$Symbol = NULL 
   
   # Drop samples having too low counts.
   tData3 = tData2[,colSums(tData2) >= cutoff]
   #tData3 = tData2
   print(paste0("number of columns of unnormalized feature table: ", ncol(tData3)))
   print(paste0("removed samples having less than ", cutoff, " reads"))
   
   x = tData3 + 1 # Deal with zeroes before submitting table to edger. Otherwise will not be able to compute norm.factor
   
   # EDGER
   y = DGEList(counts=x, remove.zeros=TRUE)
   #keep <- rowSums(cpm(y)>1) >= 2
   #y <- y[keep, , keep.lib.sizes=FALSE]
   z = edgeR::calcNormFactors(y, method="RLE")
     
   ### OLD METHOD. do not use anymore##
   #x1 = x
   #for(i in 1:ncol(x)){
   #   lib_size = z$samples[c(names(x)[i]),]$lib.size
   #   norm_fac = z$samples[c(names(x)[i]),]$norm.factors
   #   x1[,i] = round(x[,i] / (lib_size * norm_fac) * multiplier, digits=1) # downstream scripts do not support to many decimal points. So we multiply by 10000 and only keep 1 digit after .
   #}                                                 
   #print(paste0("number of columns of normalized feature table: ", ncol(x1)))
   ######################################
      
   # names
   x1 = data.frame(cpm(z), check.names=FALSE)
   x1 = round(x1, digits=1)
   x1$Symbol = row.names(x1)
   x1 = merge(x1, tTaxonomy, by.x="Symbol", by.y="Symbol", all=FALSE)
   #taxonomy_col = tData[,ncol(tData)]
   x1$taxonomy = gsub(";", "; ", x1$taxonomy)
   #x1 = cbind(x1, taxonomy_col2)
   #names(x1)[ncol(x1)] = "taxonomy"
   
   #curr_row_names = rownames(x1)
   #df = NULL
   #df = cbind(df, curr_row_names)
   #df = data.frame(df)
   #df = cbind(df, x1[,1:ncol(x1)])
   names(x1)[1] = "#FEATURE_ID"
   x2 = x1[order(rowSums(x1[,2:(ncol(x1)-1)]), decreasing=TRUE),]
   line="# normalized with edgeR RLE method"
   write(line, file=outfile ,append=FALSE)
   write.table(x2, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=TRUE)
   
   print("[DEBUG] Done Normalizing Feature table...")
} 

usage=function(errM) {
        cat("\nUsage : Rscript taxBarplotWithMappingFile.R [option] <Value>\n")
        cat("       -i        : infile\n")
        cat("       -o        : outfile\n")
        cat("       -c        : cutoff\n")
        #cat("       -m        : multiply normalized values by <n>. For shotgun metagenome, use 100000000. For 16S, typically use 10000")
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
     cutoff=ARG[i+1]
   }#else if (ARG[i] == "-m") {
    #  multiplier=ARG[i+1]
   #}
}

performNormalizationFeatures(infile, outfile, cutoff)


