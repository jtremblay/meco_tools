#!/usr/bin/env Rscript

# Function that takes a infile: (Qiime's tax summary spreadsheets) and an outdir where to
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
   
   #infile = "/home/jtrembla/Projects/testRPlots/otu_table_normalized_bacteriaArchaea_L6.txt"
   #outdir = "/home/jtrembla/Projects/testRPlots/NRCan/"
   #prefix = "test"
   #mappingFile = "/home/jtrembla/Projects/testRPlots/mapping_file.tsv"
  
   print(paste0("[DEBUG] infile: ",infile))
   print(paste0("[DEBUG] outdir: ",outdir))
   print(paste0("[DEBUG] prefix: ",prefix))
   print(paste0("[DEBUG] mapping_file: ",mappingFile))
 
   vColors = c(
      "#0000CD", "#00FF00", "#FF0000", "#808080", "#000000", "#B22222", "#DAA520", 
      "#DDA0DD", "#FF00FF", "#00FFFF", "#4682B4", "#E6E6FA", "#FF8C00", "#80008B", 
      "#8FBC8F", "#00BFFF", "#FFFF00", "#808000", "#FFCCCC", "#FFE5CC", "#FFFFCC", "#E5FFCC", 
      "#CCFFCC", "#CCFFE5", "#CCFFFF", "#CCE5FF", "#CCCCFF", "#E5CCFF", "#FFCCFF", "#FFCCE5", 
      "#FFFFFF", "#990000", "#666600", "#006666", "#330066", "#A0A0A0", "#99004C"
   )

   tData <- read.csv(file=infile, header=T, sep="\t", check.names=FALSE)
   # Hack to fix sample annoying sample renaming by Qiime
   names.orig = names(tData)
   names.mod = gsub("_", ".", names.orig)
   names.mod = gsub("-", ".", names.mod)
   names(tData) = names.mod
   ##
   tData2 <- cbind(tData, (rowSums(tData[2:ncol(tData)]))/ncol(tData) )
   tData3 <- tData2[order(-tData2[, ncol(tData2)]),]
   tData3[, ncol(tData3)] <- NULL
   tData4 <- tData3[1:20,]
   tData5 <- melt(tData4)
   order_col <- seq(1, nrow(tData5),1)
   tData5$order_col <- order_col
   tData6 <- tData5
 
   # Hack to fix sample annoying sample renaming by Qiime
   mapping_file = read.table(mappingFile, sep="\t", colClasses = "character")
   sampleID.orig = mapping_file$V1
   sampleID.mod = gsub("_", ".", sampleID.orig)
   sampleID.mod = gsub("-", ".", sampleID.mod)
   mapping_file$V1 = sampleID.mod
   ##
   header = read.table(mappingFile, sep="\t", nrows=1, comment.char="")
   colnames(mapping_file) = header
  
   for(i in 2:ncol(mapping_file)){
      cat = cbind(mapping_file[,1], mapping_file[,i])
      cat2 = data.frame(cat)
      colnames(cat2) = c("variable", names(mapping_file)[i])
      tData6 = merge(tData6, cat2, by="variable", all=FALSE)
   }

   tData7 <- tData6[with(tData6, order(order_col)),]
   tData7$order_col = NULL

   print("[DEBUG] Loaded data...")
  
   ## Sort levels of a df$variable to control the order of display on the plot.
   #tData6$Time = factor(tData6$Time, levels=sort(method="quick",levels(factor(tData6$Time))))
   #tData6$Time = factor(tData6$Time, levels=c("T0", "T5", "T42"))
  
   #generate Combinaisons (i.e. formulas...)
   #Do not generate formulas (double forumlas if mapping_file contains only 1 variable.
   print("ncol(header)")
   print(ncol(header))
   if(ncol(header) > 2){
      formulas = NULL
      for(i in 2:(ncol(header)-1)){
         for(j in (i+1):(ncol(header)-0)){
            firstVariable = as.character(header[,i])
            secondVariable = as.character(header[,j])
            formulas = c(formulas, paste(firstVariable, " ~ ", secondVariable))
         }  
      }
   }else{
      formulas = NULL
   }

   print(formulas)

   formulas_single_factor = NULL
   for(i in 2:(ncol(header))){
      currVariable = as.character(header[,i])
      formulas_single_factor = c(formulas_single_factor, paste0("~",currVariable))
   }
   
   print(formulas_single_factor)

   ## Then loop through formulas and draw plots for single factor
   for(i in formulas_single_factor){
      
      # Dynamically set font size:
      curr_variable = substring(i, 2)
      number_of_variables = length(unique(mapping_file[[curr_variable]]))
      xFontSize = -0.0125 * number_of_variables + 10.498
      yFontSize = -0.0125 * number_of_variables + 12.498
      
      if(number_of_variables > 64){
         xFontSize = -0.0125 * number_of_variables + 1.0
      }else if(number_of_variables > 32){
         xFontSize = -0.0125 * number_of_variables + 2.498
      }else if(number_of_variables > 16){
         xFontSize = -0.0125 * number_of_variables + 5.498
      }else{
         xFontSize = -0.0125 * number_of_variables + 8.498
      }
      
      prefix2 = i 
      prefix2 = gsub("~", "", prefix2)
      outfile1 <- paste0(outdir, prefix, "_", prefix2, ".pdf")
      outfile2 <- paste0(outdir, prefix, "_", prefix2, ".png")
  
      p <- ggplot(data=tData7, aes(x = variable, y=value, fill=Taxon)) + 
      facet_grid(as.formula(i), scales="free",  space = "free") + 
      geom_bar(stat="identity") + 
      theme(
         panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=1),
         axis.text.x=element_text(size=xFontSize, colour="black", angle=90, hjust=1), 
         axis.text.y=element_text(size=yFontSize, colour="black", angle=0), 
         axis.title=element_text(size=16),
         panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),
         panel.background=element_blank(),
         legend.key.size = unit(0.45, "cm"),
         legend.margin = unit(1, "cm"),
         strip.text = element_text(size=9, angle=90, hjust=0),
         strip.background =  element_blank()
      ) +  scale_fill_manual(values=vColors, breaks=rev(tData7$Taxon)) +
      labs(y="Abundance", x="Sample ID")
      pdf( file=outfile1, height=8, width=16 )
      print(p)
      dev.off()
      
      options(bitmapType='cairo') 
      png( file=outfile2, height=8, width=16, units="in", res=500)
      print(p)
      dev.off()
   }  

   ## Then loop through formulas and draw plots
   for(i in formulas){
      
      # Dynamically set font size:
      curr_variable_y = gsub(" ", "", unlist(strsplit(i, "~"))[1])
      number_of_variables_y = length(unique(mapping_file[[curr_variable_y]]))
      curr_variable_x = gsub(" ", "", unlist(strsplit(i, "~"))[2])
      number_of_variables_x = length(unique(mapping_file[[curr_variable_x]]))
      if(number_of_variables_x > 64){
         xFontSize = -0.0125 * number_of_variables_x + 1.0
      }else if(number_of_variables_x > 32){
         xFontSize = -0.0125 * number_of_variables_x + 2.498
      }else if(number_of_variables_x > 16){
         xFontSize = -0.0125 * number_of_variables_x + 5.498
      }else{
         xFontSize = -0.0125 * number_of_variables_x + 8.498
      }
      
      if(number_of_variables_y > 10){
         yFontSize = -0.0125 * number_of_variables_y + 5.498
         yFontSizeStrip = -0.0125 * number_of_variables_y + 9.498
      }else{
         yFontSize = -0.0125 * number_of_variables_y + 10.498
         yFontSizeStrip = -0.0125 * number_of_variables_y + 9.498
      }
      
      prefix2 = i
      prefix2 = gsub("  ~  ", "", prefix2)
      outfile1 <- paste0(outdir, prefix, "_", prefix2, ".pdf")
      outfile2 <- paste0(outdir, prefix, "_", prefix2, ".png")
    
      p <- ggplot(data=tData7, aes(x = variable, y=value, fill=Taxon)) + 
         facet_grid(as.formula(i), scales="free",  space = "free_x") + 
         geom_bar(stat="identity") + 
         theme(
            panel.border=element_rect(fill=NA, linetype="solid", colour = "black", size=1),
            axis.text.x=element_text(size=xFontSize, colour="black", angle=90, hjust=1), 
            axis.text.y=element_text(size=yFontSize, colour="black", angle=0), 
            axis.title=element_text(size=16),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.background=element_blank(),
            legend.key.size = unit(0.45, "cm"),
            legend.margin = unit(1, "cm"),
            strip.text.x = element_text(size=xFontSize, angle=90, hjust=0),
            strip.text.y = element_text(size=yFontSizeStrip, angle=0, hjust=0),
            strip.background =  element_blank()
      ) +  scale_fill_manual(values=vColors, breaks=rev(tData7$Taxon)) +
      labs(y="Abundance", x="Sample ID")
      pdf( file=outfile1, height=8, width=16 )
      print(p)
      dev.off()

      png( file=outfile2, height=8, width=16, units="in", res=500)
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

