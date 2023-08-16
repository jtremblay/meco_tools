#!/usr/bin/env Rscript

# Function that generates pcoas based on abundance matrix. Typically for MG and MT expr. data.
# To write plot files. julien.tremblay@nrc-cnrc.gc.ca
generatePcoAPlots <- function(indir, mappingFile, outdir) {
  library(ggplot2)
  library(scales)
  library(reshape2)
  library(grid)
  library(gridExtra)
  library(tools)
  library(labdsv)
  library(data.table)
   
  options(stringsAsFactors = FALSE)
   
  #root        = "~/Projects/test_MT"
  #infile      = paste0(root, "/merged_gene_abundance_cpm.tsv")
  #mappingFile = paste0(root, "/mapping_file.tsv")
  #outdir      = paste0(root, "/out_pcoa/")
  
  dir.create(file.path(outdir), showWarnings=FALSE)
  
  print(paste0("[DEBUG] infile: ", infile))
  print(paste0("[DEBUG] outdir: ", outdir))
  print(paste0("[DEBUG] mapping file: ", mappingFile))

  vColors = c(
     "#0000CD", "#00FF00", "#FF0000", "#808080", "#000000", "#B22222", 
     "#40E0D0", "#DAA520", "#DDA0DD", "#FF00FF", "#00FFFF", "#4682B4", 
     "#008000", "#E6E6FA", "#FF8C00", "#80008B", "#8FBC8F", "#00BFFF", 
     "#FFFF00", "#808000", "#FFCCCC", "#FFE5CC", "#FFFFCC", "#E5FFCC", 
     "#CCFFCC", "#CCFFE5", "#CCFFFF", "#CCE5FF", "#CCCCFF", "#E5CCFF", 
     "#FFCCFF", "#FFCCE5", "#FFFFFF", "#990000", "#666600", "#006666",
     "#330066","#A0A0A0","#99004C"
  )

  tData = data.frame(fread(infile, header=TRUE, sep="\t"), check.names=FALSE, row.names=1)
  print(head(tData))

  dist_obs = dist(t(tData))
  dist_matrix = as.matrix(dist_obs)
  coords = pco(dist_obs, k=ncol(dist_matrix)-1)
  coord_matrix = coords$points
  coord_df = data.frame(coord_matrix)
  eig = coords[[2]]
  options("scipen"=100, "digits"=4)
  percent_values = as.numeric((eig / sum(eig))*100)
  percent_values = round(percent_values, digits=2)
  percent1 = percent_values[1]
  percent2 = percent_values[2]
  tData2 <- coord_df[,1:3]
  tData3 = data.frame(row.names(tData2))
  tData4 = cbind(tData3, tData2)
  colnames(tData4) = c("variable", "D1", "D2", "D3")
       
  mapping_file = read.table(mappingFile, sep="\t")
  header = read.table(mappingFile, sep="\t", nrows=1, comment.char="")
  colnames(mapping_file) = header
  
  for(i in 2:ncol(mapping_file)){
    cat = cbind(mapping_file[,1], mapping_file[,i])
    cat2 = data.frame(cat)
    colnames(cat2) = c("variable", names(mapping_file)[i])
    tData4 = merge(tData4, cat2, by="variable", all=FALSE)
  }

  print("[DEBUG] Loaded data...")

  ## Sort levels of a df$variable to control the order of display on the plot.
  #generate Combinaisons
  tFirstVariable = NULL
  tSecondVariable = NULL
  if(ncol(header) == 2){
     firstVariable = as.character(header[,2])
     secondVariable = as.character(header[,2])
     tFirstVariable = c(tFirstVariable, paste(firstVariable))
     tSecondVariable = c(tSecondVariable, paste(secondVariable))
  }else{
    for(i in 2:(ncol(header)-1)){
      for(j in (i+1):(ncol(header)-0)){
        firstVariable = as.character(header[,i])
        secondVariable = as.character(header[,j])
        tFirstVariable = c(tFirstVariable, paste(firstVariable))
        tSecondVariable = c(tSecondVariable, paste(secondVariable))
       } 
     }
  }

  .e <- environment()
  plots = list()
  for(i in 1:length(tFirstVariable)){
    prefix2 = paste0(tFirstVariable[i], "_", tSecondVariable[i])
    first = tFirstVariable[i]
    second = tSecondVariable[i]   
   
    p <- ggplot(environment=.e, data=tData4, aes(x=D1, y=D2, color=as.character(tData4[[tFirstVariable[i]]]), shape=as.character(tData4[[tSecondVariable[i]]]))) + 
       geom_point(size=5) + geom_point(colour="grey90", size = 1.5) +
       theme(
          panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=1.4),
          axis.text.x=element_text(size=15, colour="black"), 
          axis.text.y=element_text(size=15, colour="black"), 
          axis.title=element_text(size=16),
          axis.ticks.length=unit(0.2,"cm"),
          axis.ticks = element_line(colour = 'black', size = 0.5),
          panel.grid.minor=element_blank(),
          panel.background=element_blank(),
          legend.key.size = unit(0.45, "cm"),
          legend.margin = unit(1, "cm")
       ) + scale_colour_manual(values=vColors, name = first) + 
       geom_hline(aes(yintercept=0), size=0.2, linetype="dotted") + 
       xlab(paste0("PCo1 (", percent1,"%)")) + 
       ylab(paste0("PCo2 (", percent2,"%)")) +
           #scale_colour_discrete(name = first) +
       scale_shape_discrete(name = second)
       #geom_point(size=5) + geom_point(colour="grey90", size = 1.5) + geom_text(aes(label=variable), vjust=1.5) +  
        
    y = paste0(outdir, "/")
    outfile_base = file_path_sans_ext(y)
    outfile1 = paste0(outfile_base, prefix2, "_pcoa.pdf")
    outfile2 = paste0(outfile_base, prefix2, "_pcoa.jpeg")
    pdf( file=outfile1, height=6, width=10 )
    #do.call(grid.arrange, c(plots, ncol=2))
    print(p)
    dev.off()
    jpeg( file=outfile2, height=6, width=10, units="in", res=500)
    #do.call(grid.arrange, c(plots, ncol=2))
    print(p)
    dev.off()
    
    prefix3 = paste0(tSecondVariable[i], "_", tFirstVariable[i])
    p <- ggplot(environment=.e, data=tData4, aes(x=D1, y=D2, color=as.character(tData4[[tSecondVariable[i]]]), shape=as.character(tData4[[tFirstVariable[i]]]))) + 
       geom_point(size=5) + geom_point(colour="grey90", size = 1.5) +
       theme(
          panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=1.4),
          axis.text.x=element_text(size=15, colour="black"), 
          axis.text.y=element_text(size=15, colour="black"), 
          axis.title=element_text(size=16),
          axis.ticks.length=unit(0.2,"cm"),
          axis.ticks = element_line(colour = 'black', size = 0.5),
          panel.grid.minor=element_blank(),
          panel.background=element_blank(),
          legend.key.size = unit(0.45, "cm"),
          legend.margin = unit(1, "cm")
       ) + scale_colour_manual(values=vColors, name=second) + 
       geom_hline(aes(yintercept=0), size=0.2, linetype="dotted") + 
       xlab(paste0("PCo1 (", percent1,"%)")) + 
       ylab(paste0("PCo2 (", percent2,"%)")) +
       #scale_colour_discrete(name = first) +
       scale_shape_discrete(name = first)
    #geom_point(size=5) + geom_point(colour="grey90", size = 1.5) + geom_text(aes(label=variable), vjust=1.5) +  
    
    #print(p)
    y = paste0(outdir, "/")
    outfile_base = file_path_sans_ext(y)
    outfile1 = paste0(outfile_base, prefix3, "_pcoa.pdf")
    outfile2 = paste0(outfile_base, prefix3, "_pcoa.jpeg")
    pdf( file=outfile1, height=6, width=10 )
    #do.call(grid.arrange, c(plots, ncol=2))
    print(p)
    dev.off()
    jpeg( file=outfile2, height=6, width=10, units="in", res=500)
    #do.call(grid.arrange, c(plots, ncol=2))
    print(p)
    dev.off()
    
  }
  print("[DEBUG] Printed plots...")
  
}

usage=function(errM) {
        cat("\nUsage : Rscript plotAbundancePcoa.R [option] <Value>\n")
        cat("       -i        : indir\n")
        cat("       -p        : prefix\n")
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
	}else if (ARG[i] == "-m") {
		mappingFile=ARG[i+1]
	}else if (ARG[i] == "-o") {
		outdir=ARG[i+1]
	}
}

generatePcoAPlots(infile, mappingFile, outdir)
