#!/usr/bin/env Rscript

# Function that takes a Feature table and generate a heatmap.
# National Research Council Canada
# Genomics and Microbiomes
# Author: julien.tremblay@nrc-cnrc.gc.ca
# To write plot files.
makeplots <- function(Feature_table, outdir, prefix, n, mappingFile) {

  options(stringsAsFactors = FALSE)
  library(data.table)
  library(pheatmap)

  filename = Feature_table
  outfilePng = paste0(outdir,'/', prefix, '.png')
  outfilePdf = paste0(outdir,'/', prefix, '.pdf')

  #======
  #filename = "~/Projects/heatmap_R/feature_table_filtered_rarefied.tsv"
  #mappingFile = "~/Projects/heatmap_R/mapping_file.txt"
  #n=15
  #======
  
  # Process mapping file
  mapping_file = read.table(mappingFile, sep="\t")
  header = read.table(mappingFile, sep="\t", nrows=1, comment.char="")
  header = gsub("\\#", "", header)
  colnames(mapping_file) = header
  d = sapply(mapping_file, as.character)
  d1 = data.frame(d)
  rownames(d1) = d1$SampleID
  d1$SampleID = NULL
  mapping_file = d1
      
  # Process feature table
  data = data.frame(fread(filename, sep="\t", header=FALSE, skip="#FEATURE_ID", colClasses="character"), check.names=FALSE)
  data = as.matrix(data)
  data2=data[2:nrow(data),3:ncol(data)-1]
  colLabels = data[1,]
  rowLabels = data[,1]
  #print(colLabels)
  lineages = data[,ncol(data)]
  data3 = data.frame(data2,  stringsAsFactors = FALSE)
  names(data3) = colLabels[3:length(colLabels)-1]
  IDAndLineages = paste(rowLabels, lineages)
  rownames(data3) = IDAndLineages[2:length(IDAndLineages)]
 
  # Continue processing mapping_file 
  #row.names(mapping_file) = mapping_file$SampleID
  #mapping_file$SampleID = NULL
  #print(mapping_file)   
 
  mat = data.matrix(data3)
  matOrdered = mat[rev(order(rowSums(mat))),]


  if(n == 0){
    matOrderedSubset = matOrdered
  }else{
    curr_n = nrow(matOrdered)
    n = as.numeric(n)
    #print(n) #40
    #print(curr_n) #54
    if(curr_n < n){
      quit()
    }
    matOrderedSubset = matOrdered[1:n,]
  }
  matOrderedSubset = matOrderedSubset + 1
  matOrderedSubset = log2(matOrderedSubset)
  
  #====
  #pheatmap(matOrderedSubset, fontsize_row=6, fontsize_col=6, annotation=mapping_file, clustering_method="average")
  #====
  
  #pheatmap(matOrderedSubset, file=outfilePng, fontsize_row=6, fontsize_col=6,cellwidth=10, cellheight=6, annotation=mapping_file, clustering_method="average")
  #pheatmap(matOrderedSubset, file=outfilePdf, fontsize_row=6, fontsize_col=6,cellwidth=10, cellheight=6, annotation=mapping_file, clustering_method="average")  
  font_size = 6
  font_size_col = 6
  if(n == 20){
    font_size = 11
    font_size_col = 9
    w = 24
  } else if(n == 40){
    font_size = 9
    font_size_col = 8
    w = 20
  } else if(n == 60){
    font_size = 5
    w = 16
  } else if(n == 80){
    font_size = 3
    w = 12
  } else if(n == 100){
    font_size = 2
    w = 10
  }
   
  # Then adjust according to number of samples...
  n_samples = ncol(matOrderedSubset)
  w = (n_samples *13 ) / 30
  w_cell = (n_samples * 6) / 30

  # For very low sample number
  if(n_samples < 10){
    w = (n_samples *13 ) / 1.7
    w_cell = (n_samples * 6) / 1
  }
  
   # For high sample number
  if(n_samples > 200){
    w = (n_samples *13 ) / 170
    w_cell = (n_samples * 6) / 1000 
    font_size_col = 0.5
  }

  print(n_samples)
  print(w)
  print(w_cell)
  #print(mapping_file)
  options(bitmapType='cairo')
  pheatmap(main="Feature heatmap (log2)", matOrderedSubset, file=outfilePng, fontsize_row=font_size, fontsize_col=font_size_col, cellwidth=(w_cell),  width=w,  annotation=mapping_file, clustering_method="average")
  pheatmap(main="Feature heatmap (log2)", matOrderedSubset, file=outfilePdf, fontsize_row=font_size, fontsize_col=font_size_col,  cellwidth=(w_cell),  width=w,  annotation=mapping_file, clustering_method="average")  
} 

usage=function(errM) {
        cat("\nUsage : Rscript Featureheatmap.R [option] <Value>\n")
        cat("       -t        : Feature table\n")
        cat("       -o        : outdir\n")
        cat("       -p        : prefix\n")
        cat("       -n        : Number of rows to include in heatmap\n")
        cat("       -m        : Mapping file\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 5) {
   usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
   if(ARG[i] == "-t") {
      Feature_table=ARG[i+1]
   }else if (ARG[i] == "-o") {
      outdir=ARG[i+1]
   }else if (ARG[i] == "-p") {
      prefix=ARG[i+1]
   }else if (ARG[i] == "-n") {
     n=ARG[i+1]
   }else if (ARG[i] == "-m") {
     mappingFile=ARG[i+1]
   }
}

makeplots(Feature_table, outdir, prefix, n, mappingFile)
