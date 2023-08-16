#!/usr/bin/env Rscript

# Function that takes multiples gene count bed files and 
# merges them into on file (matrix). Intended for the mgs_augmented_assembly pipeline
# Author: Julien Tremblay - jtremblay514@gmail.com

options(stringsAsFactors = FALSE)
mergeTables <- function(infiles, outfileRaw, outfileCpm, type, skipNorm) {

  library(edgeR)
  library(data.table)

  if(is.null(skipNorm)){
     skipNorm = "FALSE"
  }
  if(skipNorm == "FALSE"){
    skipNorm = FALSE
  }else{
    skipNorm = TRUE
  }

  print("skipNorm:")
  print(skipNorm)

  index = NULL
  symbol_index = NULL
  if(type == "contigs"){
    index = 4 
    symbol_index = 1
  }else if(type == "genes"){
    index = 5
    symbol_index = 4
  }else{
    stop("-t arg must be either 'contigs' or 'genes'")
  }
  
  files = unlist(strsplit(infiles, ","))
  
  df = NULL
  #dfNormalized = NULL
  
  #tData = read.table(files[1], sep="\t", header=FALSE, comment.char="")
  tData = data.frame(fread(files[1], sep="\t", header=FALSE, showProgress=FALSE), check.names=FALSE)

  symbol = as.character(tData[,symbol_index])
  #geneName = as.character(tData[,2])
 
  df = cbind(df, symbol)
  #dfNormalized = cbind(dfNormalized, symbol)
  
  #print("df:")
  #print(head(df))

  if(type == "contigs"){
    colNames = c("contig_id")
  }else if(type == "genes"){
    colNames= c("gene_id")
  }
  
  for(i in 1:length(files)){
    print(paste0("Reading file ", files[i]))
    currFile = files[i];
    currName = basename(currFile)
    currName = gsub(".cov", "", currName)
    colNames = append(colNames, currName)
    #currTable = read.table(currFile, sep="\t", header=FALSE, comment.char="")
    currTable = data.frame(fread(files[i], sep="\t", header=FALSE, showProgress=FALSE), check.names=FALSE)
    #print("currTable:")
    #print(head(currTable))
    #stop()
    #length = (currTable[,3] - currTable[,2]) / 1000 #To get gene lengths in kb.
    #totalReads = sum(currTable[,index])
    #normalizedCount = (currTable[,4] / (length * 1000)) # report mapped reads per kb.
    rawCount = (currTable[,index]) # report mapped reads per kb.
    ###RPKM = (currTable[,4] / (length * 1000) / (totalReads / 1000000)) # report mapped reads per kb per million reads mapped (RPKM).
    #RPKM = (1000000 * currTable[,index]) / (length * totalReads) # report mapped reads per kb per million reads mapped (RPKM).
 
    df = cbind(df, rawCount)
    #dfNormalized = cbind(dfNormalized, RPKM)
    #print(paste0("Total reads:", totalReads))
  }
  df = df[!duplicated(df), ]
  df = data.frame(df)
  #dfNormalized = data.frame(dfNormalized)
  names(df) = colNames
  #names(dfNormalized) = colNames
  df[is.na(df)] = 0
  #dfNormalized[is.na(dfNormalized)] = 0
  
  write.table(df, file=outfileRaw, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  #write.table(dfNormalized, file=outfileNormalized, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
 
  if(isFALSE(skipNorm)){
      # Then generate cpms
      row.names(df) = df[,1]
      df$cluster = NULL
      df$contig_id = NULL
      df$gene_id = NULL
      tmpGroup = c(rep("A", ncol(df)))
      
      for(i in 1:ncol(df)){
         df[,i] = as.numeric(df[,i])
      }

      #row.names(tData) = tData$V1
      #tData$V1 = NULL
      #print(head(tData))

      ## edgeR. Load table, remove low cpm. 
      #tData = tData[which(rowSums(tData) > 0),]
      y <- DGEList(df, remove.zeros=FALSE)
      y <- calcNormFactors(y, method="TMM")
      cpms = cpm(y)
      cpms = round(cpms, digits=3)
      #keep <- rowSums(cpm(y)>1) >= 2
      #y <- y[keep, , keep.lib.sizes=FALSE]
      #cpms = apply(df,2, function(x) (round( ((x/(sum(x)+1))*1000000), digits=3)))
      #cpms = round(cpm(y1), digits=3)
      
      
      #cpms = round(cpms, digits=3)
      #write.table(cpms, outfileCpm, quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)
      fwrite(as.data.frame(cpms), outfileCpm, row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)
      
      print("Done merging and normalizing tables...")
    } 
} 

usage=function(errM) {
  cat("\nUsage : Rscript mergeRawCounts.R [option] <Value>\n")
  cat("       -i        : List of files separated by a comma\n")
  cat("       -o        : outfile raw counts\n")
  cat("       -c        : outfile CPMs\n")
  cat("       -t        : 'genes' or 'contigs'\n")
  cat("       -s        : skip normalization\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 5) {
  usage("missing arguments")
}

skipNorm = NULL
## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-i") {
    infiles = ARG[i+1]
  } else if (ARG[i] == "-o") {
    outfile = ARG[i+1]
  #} else if (ARG[i] == "-n") {
  #  outfileNormalized = ARG[i+1]
  } else if (ARG[i] == "-c") {
    outfileCpm = ARG[i+1]
  } else if (ARG[i] == "-t") {
    type = ARG[i+1]
  } else if (ARG[i] == "-s") {
    skipNorm = ARG[i+1]
  }
}

mergeTables(infiles, outfile, outfileCpm, type, skipNorm)
