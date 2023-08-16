#!/usr/bin/env Rscript

# Function that takes a infile: (Qiime's PCoA coordinates file) and an outdir where to
# National Research Council Canada
# Genomics and Microbiomes
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
# To write plot files.
generatePlots <- function(infile, outdir, prefix, mappingFile) {
   library(ggplot2)
   library(scales)
   library(reshape2)
   library(grid)
   options(stringsAsFactors = FALSE)
   
   #infile = "~/Projects/testRPlots/unweighted_unifrac_feature_table_final_normalized_coords.tsv"
   #outdir = "~/Projects/testRPlots/"
   #prefix = "pca-test"
   #mappingFile = "~/Projects/testRPlots/mapping_file_DFO.tsv"
  
   print(paste0("[DEBUG] infile: ",infile))
   print(paste0("[DEBUG] outdir: ",outdir))
   print(paste0("[DEBUG] prefix: ",prefix))
   print(paste0("[DEBUG] prefix: ",mappingFile))

  
  #vColors = c(
  #  '#990099',
  #  '#009900',
  #  '#0099FF',
  #  '#FF6600',
  #  '#000000',
  #  '#999999',
  #  '#9999FF',
  #  '#FF0000'
  #)
   vColors = c(
     '#00FF00',
     '#FF8080',
     '#FF00FF',
     '#0000FF',
     '#00CCFF',
     '#CCFFFF',
     '#CCFFCC',
     '#99CCFF',
     '#CC99FF',
     '#FFCC99',
     '#3366FF',
     '#33CCCC',
     '#99CC00',
     '#FF99CC',
     '#FFCC00'
   )
  
   # Directly read into table before anything else.
   coords_file = file(infile)
   open(coords_file)
   percentVar = read.table(coords_file, skip=4, nrow=1)
   percent1 = percentVar[[1]] * 100
   percent2 = percentVar[[2]] * 100
   percent1 = formatC(round(percent1, 2), big.mark=",",drop0trailing=TRUE, format="f")
   percent2 = formatC(round(percent2, 2), big.mark=",",drop0trailing=TRUE, format="f")
   close(coords_file)
   
   # Then process taxonomy
   coords_file = file(infile)
   open(coords_file)
   # Read starting from beginning of coordinates, fill last rows with NAs
   tData = read.table(coords_file, skip=9, fill=TRUE)
   # Then remove NAs
   tData = na.omit(tData)
   # Keep 3 first rows only.
   tData2 = tData[,1:4]
   for(i in 2:ncol(tData2)){
      tData2[,i] = as.numeric(tData2[,i]) 
   }
   tData3 = tData2
   #row.names(tData) = tData[,1]
   #tData$V1 = NULL
   colnames(tData3) = c("variable", "D1", "D2", "D3")
   close(coords_file)
 
   mapping_file = read.table(mappingFile, sep="\t")
   header = read.table(mappingFile, sep="\t", nrows=1, comment.char="")
   colnames(mapping_file) = header
  
   for(i in 2:ncol(mapping_file)){
      cat = cbind(mapping_file[,1], mapping_file[,i])
      cat2 = data.frame(cat)
      colnames(cat2) = c("variable", names(mapping_file)[i])
      tData3 = merge(tData3, cat2, by="variable", all=FALSE)
   }
  
   # Make sure everyting is in int.

   print("[DEBUG] Loaded data...")

   ## Sort levels of a df$variable to control the order of display on the plot.
   #tData6$Time = factor(tData6$Time, levels=sort(method="quick",levels(factor(tData6$Time))))
   #tData6$Time = factor(tData6$Time, levels=c("T0", "T5", "T42"))
  
   #generate Combinaisons
   #header$V2 = NULL
   tFirstVariable = NULL
   tSecondVariable = NULL
   if(ncol(header) > 2){
      for(i in 2:(ncol(header)-1)){
         for(j in (i+1):(ncol(header)-0)){
            firstVariable = as.character(header[,i])
            secondVariable = as.character(header[,j])
            tFirstVariable = c(tFirstVariable, paste(firstVariable))
            tSecondVariable = c(tSecondVariable, paste(secondVariable))
         } 
      }
   }else{
      tFirstVariable = c(as.character(header[,2]))
      tSecondVariable = tFirstVariable
   }
  
  # Function to manage number of tick for x-y axes
  #number_ticks <- function(n) {function(limits) pretty(limits, n)}

   .e <- environment()
   for(i in 1:length(tFirstVariable)){
      prefix2 = paste0(tFirstVariable[i], "_", tSecondVariable[i])
      first = tFirstVariable[i]
      second = tSecondVariable[i]
      outfile1 <- paste0(outdir, prefix, "_", prefix2, ".pdf")
      outfile2 <- paste0(outdir, prefix, "_", prefix2, ".png")
    
      p <- ggplot(environment=.e, data=tData3, aes(x=D1, y=D2, color=as.character(tData3[[first]]), shape=as.character(tData3[[second]]))) + 
      geom_point(size=5) +  geom_point(colour="grey90", size = 1.5) + 
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
      ) + scale_fill_manual(values=vColors) + 
      geom_hline(aes(yintercept=0), size=0.2) + 
      xlab(paste0("PCo1 (", percent1,"%)")) + 
      ylab(paste0("PCo2 (", percent2,"%)")) +
      scale_colour_discrete(name = first) +
      scale_shape_discrete(name = second) 
    
     pdf( file=outfile1, height=6, width=8 )
     print(p)
     dev.off()
    
     options(bitmapType='cairo')
     png( file=outfile2, height=6, width=8, units="in", res=300)
     print(p)
     dev.off()
  }

	
	print("[DEBUG] Printed plots...")
} 

usage=function(errM) {
        cat("\nUsage : Rscript taxBarplotWithMappingFile.R [option] <Value>\n")
        cat("       -i        : infile\n")
        cat("       -o        : outdir\n")
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
		infile=ARG[i+1]
	} else if (ARG[i] == "-o") {
		outdir=ARG[i+1]
	}else if (ARG[i] == "-p") {
		prefix=ARG[i+1]
	}else if (ARG[i] == "-m") {
		mappingFile=ARG[i+1]
	}
}

generatePlots(infile, outdir, prefix, mappingFile)

