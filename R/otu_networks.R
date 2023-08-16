#!/usr/bin/env Rscript

# Function that generate networks plot.
generatePlots <- function(infile, otu_table, outdir) {
   library(ggplot2)
   library(scales)
   library(reshape2)
   library(grid)
   library(igraph)
   library(scales)
   library(data.table)
   library(RColorBrewer)
   options(stringsAsFactors = FALSE)
   
   #infile = "~/Projects/testRPlots/unweighted_unifrac_otu_table_final_normalized_coords.tsv"
   infile = "~/Projects/networks/otu_test_cor_sparcc.tsv"
   otu_table = data.frame(fread("~/Projects/networks/otu_table_filtered_bacteriaArchaea_normalized.tsv"), check.names=FALSE)
   outdir = "/home/jtrembla/Projects/networks/"
   low_threshold = -0.63
   high_threshold = 0.63
   
   #prefix = "network-test"
   #mappingFile = "~/Projects/testRPlots/mapping_file_DFO.tsv"
  
   print(paste0("[DEBUG] infile: ",infile))
   print(paste0("[DEBUG] outdir: ",outdir))

   outfile1 <- paste0(outdir, prefix, "_", ".pdf")
   outfile2 <- paste0(outdir, prefix, "_", ".jpeg")

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
   
   # chomp and chop function
   trim <- function (x) gsub("^\\s+|\\s+$", "", x)
   trim2 <- function (x) gsub("; ", ";", x)
   trim3 <- function (x) gsub(";$", "", x)
   
   # Function for padding
   len=6
   na.pad <- function(x,len){x[1:len]}
   makePaddedDataFrame <- function(l,...){
      maxlen <- max(sapply(l,length))
      data.frame(lapply(l,na.pad,len=maxlen),...)
   }
   
   # function for placing vertice label outside of ploting area
   radian.rescale <- function(x, start=0, direction=1) {
      c.rotate <- function(x) (x + start) %% (2 * pi) * direction
      c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
   }
   

   df = data.frame(fread(infile, header=TRUE), check.names=FALSE)
   row.names(df) = df[,1]
   df[,1] = NULL
   cor.mat = as.matrix(df)
   cor.mat <- ifelse(cor.mat > low_threshold & cor.mat < high_threshold, 0, cor.mat) # keep values > abs(0.3)
   cor.mat <- ifelse(cor.mat == 1, 0, cor.mat) # keep values > abs(0.3)1 # remove same OTU comparison.
   cor.mat = as.data.frame(as.table(cor.mat))
   cor.mat = cor.mat[cor.mat[,3]!=0,]
      
   ## Then open otu table to associate lineage with each otu #.
   
   otu_table2 = data.frame(otu_table[,1], otu_table[,ncol(otu_table)])
   names(otu_table2) = c("OTU_ID", "lineage")
   # Add otu sums in fraction 
   sums = rowSums(otu_table[,2:(ncol(otu_table)-1)])
   otu_table2$size = sums
   otu_table2$perc = (otu_table2$size / sum(otu_table2$size)) * 100
   otu_table2$size = NULL
   
   # create a size table for later merge.
   size_table = NULL
   size_table = data.frame(cbind(otu_table2$OTU_ID, otu_table2$perc))
   names(size_table) = c("OTU_ID", "perc")
   
   otu_table2$lineage = trim(otu_table2$lineage)
   otu_table2$lineage = trim2(otu_table2$lineage)
   otu_table2$lineage = trim3(otu_table2$lineage)
   tmp = strsplit(as.character(otu_table2$lineage), ";", fixed=TRUE)
   tmp2 = makePaddedDataFrame(tmp)
   tmp3 = data.frame(t(tmp2), row.names=NULL)
   names(tmp3) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
   otu_table2 = cbind(otu_table2, tmp3)
   # Here choose tax level
   otu_table3 = data.frame(otu_table2$OTU_ID, otu_table2$Genus)
   names(otu_table3) = c("OTU_ID", "Order")
   otu_table3 = merge(otu_table3, size_table, by="OTU_ID")
   otu_table3$perc.x = NULL
   names(otu_table3) = c("OTU_ID", "Order", "perc")
   
   # Then process and plot network.
   g = graph_from_data_frame(cor.mat, directed=TRUE)
   V(g)$label <- V(g)$name
   V(g)$degree <- degree(g)
   #V(g)$size <- log2(degree(g)*10)
   #V(g)$size <- degree(g)/5
   V(g)$size <- otu_table3[otu_table3$OTU_ID %in% curr_otus, 3] * 2
   V(g)$size[V(g)$size < 2] = 2
   #V(g)$size
   neg_pal = colorRampPalette(c('orangered','darkred'))
   pos_pal = colorRampPalette(c('lightgreen','darkgreen'))
   cols_neg <- neg_pal(10)[as.numeric(cut( (E(g)[Freq<0]), breaks=10))]
   cols_pos <- pos_pal(10)[as.numeric(cut( (E(g)[Freq>0]), breaks=10))]
   E(g)[Freq < 0 ]$color = cols_neg
   E(g)[Freq > 0 ]$color = cols_pos
   E(g)$weight = 1.5
   
   # Color data frame.
   curr_otus = V(g)$label
   df_colors = otu_table3[otu_table3$OTU_ID %in% curr_otus,]
   unique_orders = unique(df_colors$Order)
   df_colors_unique = data.frame(unique_orders)
   df_colors_unique = cbind(df_colors_unique, vColors[1:nrow(df_colors_unique)])
   names(df_colors_unique) = c("Order", "color")
   df_colors2 = merge(df_colors, df_colors_unique, by="Order")
   df_colors2 = df_colors2[match(curr_otus, df_colors2$OTU_ID),]
   V(g)$color = alpha(df_colors2$color, 1.0)
      
   par(mai=c(1,1,1,1))
   lab.locs <- radian.rescale(x=1:nrow(df_colors2), direction=-1, start=0)
   layout = layout.circle(g)
   plot(g, layout=layout, edge.arrow.size=0, edge.width=E(g)$weight, vertex.label.dist=1, vertex.label.degree=lab.locs)
   layout = layout.kamada.kawai(g)
   plot(g, layout=layout, edge.arrow.size=0.01, edge.width=E(g)$weight, edge.curved=FALSE, vertex.label.dist=0.4)
   legend('topright', legend=df_colors_unique$Order, col=df_colors_unique$color, cex=1, pt.cex=3, pch=20)

   #layout.circle(graph, params)
   #layout.sphere(graph, params)
   #layout.fruchterman.reingold(graph, ..., dim=2, params)
   #layout.kamada.kawai(graph, ..., dim=2, params)
   #layout.spring(graph, ..., params)
   #layout.reingold.tilford(graph, ..., params) # straight lines.
   #layout.fruchterman.reingold.grid(graph, ..., params)
	
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

