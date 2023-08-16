#!/usr/bin/env Rscript

# Function that will draw heatmap plot of supplied gene abundance table.
# Julien Tremblay => julien.tremblay@nrc-cnrc.gc.ca
drawHeatmap <- function(indir, mappingFile, prefix) {
   #library(heatmap3)
   library(RColorBrewer)
   library(pheatmap)
   library(edgeR)
   library(tools)
   
   options(stringsAsFactors = FALSE)
   
   #infile = "~/Projects/metagenome_stats/DDA_kegg_parsed.tsv"
   #outdir = "~/Projects/metagenome_stats/"
   #mappingFile = "~/Projects/metagenome_stats/mapping_file2.tsv"
   
   # Mapping file
   mapping_file = read.table(mappingFile, sep="\t", colClasses = "character")
   header = read.table(mappingFile, sep="\t", nrows=1, comment.char="")
   colnames(mapping_file) = header
   
   files = list.files(path=indir, pattern=paste0(prefix, "$"), full.names=TRUE, recursive=TRUE)
   print(files)
   for(y in files){
      print(y)
      outfile = file_path_sans_ext(y)
      outfile = paste0(outfile, "_heatmap.pdf")
      print(outfile)
      # Read abundance table
      tData <- read.csv(file=y, header=T, skip=0, sep="\t", check.names=FALSE)
      row.names(tData) = tData[,1] 
      tData$cluster = NULL
      col_names = colnames(tData[,1:ncol(tData)])
      mapping_file_parsed = mapping_file[mapping_file[,1] %in% col_names,]
      
      df = NULL
      col_names = colnames(tData[,2:ncol(tData)])
      df = cbind(df, col_names)
      df = data.frame(df)
      #df = as.character(df)
      vec = mapping_file_parsed[df[,1] %in% mapping_file[,1],]
      vec$CombinedFactors1 = paste(vec[,2:(ncol(mapping_file))])
      vec = vec$CombinedFactors1
      print(levels(factor(vec)))
 
      #df_col = NULL
      #df_col = levels(factor(vec))
      #df_col = data.frame(df_col)
      #myColors = brewer.pal(length(levels(factor(vec))),  "RdBu")
      #df_col = cbind(df_col, myColors)
      #myColors = df_col$myColors[match(vec, df_col$df_col)]
      if(nrow(tData) < 2){
         next;
      }
      # Normalize edger
      #x = tData + 1 
      x1 = tData
 
      # EDGER
      #y = edgeR::DGEList(counts = x, remove.zeros = TRUE)
      #z = edgeR::calcNormFactors(y, method="RLE")
      
      #x1 = x 
      #for(i in 1:ncol(x)){
      #  lib_size = z$samples[c(names(x)[i]),]$lib.size
      #  norm_fac = z$samples[c(names(x)[i]),]$norm.factors
      #  x1[,i] = x[,i] / (lib_size * norm_fac)
      #} 
      #print("Done normalizing table with Edger")  
      
      # Arrange final table
      #factors = mapping_file_parsed[df[,1] %in% mapping_file_parsed[,1],]
      #factors = factors$CombinedFactors1
      factors = vec
      names = mapping_file_parsed[df[,1] %in% mapping_file_parsed[,1], 1]
      df_names = NULL
      df_names = cbind(df_names, names)
      df_names = data.frame(df_names)
      df_names = cbind(df_names, factors)
      
      curr_col_names = df_names[match(colnames(x1), df_names$names),1]
      curr_factor_names = df_names[match(curr_col_names, df_names$names),2]
      
      combined_names = paste0(curr_col_names, "<=>", curr_factor_names)
      #colnames(x) = combined_names
      colnames(x1) = combined_names

      # Convert to cpm
      x2 = x1 * 1000000
      x2 = x2 + 0.0001 # Deal with zeros before transforming to log scale
      #x2 = log10(x2[1:100,])
      x2 = log2(x2)
      
      pheatmap(x2, 
         filename=outfile, 
         #colSideColors=myColors, 
         fontsize_row=2,
         fontsize_col=3,
         cellwidth=5,
         cellheight=2,
         clustering_method="complete"
      )
   }
}

usage=function(errM) {
        cat("\nUsage : Rscript heatmap.R [option] <Value>\n")
        cat("       -i        : indir\n")
        #cat("       -o        : outdir\n")
        cat("       -p        : prefix\n")
        cat("       -m        : mapping_file\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 4) {
    usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-i") {
        indir=ARG[i+1]
    #} else if (ARG[i] == "-o") {
       # outdir=ARG[i+1]
    }else if (ARG[i] == "-m") {
        mappingFile=ARG[i+1]
    }else if (ARG[i] == "-p") {
      prefix=ARG[i+1]
    }
}

drawHeatmap(indir, mappingFile, prefix)
