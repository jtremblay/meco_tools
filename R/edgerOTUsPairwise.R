#!/usr/bin/env Rscript

#  Julien Tremblay => julien.tremblay@nrc-cnrc.gc.ca
performEdgerPairwiseOTUs <- function(infile, outdir, designFile, pvalue, fdr, logfc) {
   library(edgeR)
   library(tools)
   library(data.table)

   options(stringsAsFactors = FALSE)
   
   #root       = "/home/jtrembla/Projects/test_edger_OTUs"
   #infile     = paste0(root, "/otu_table_filtered_bacteriaArchaea.tsv")
   #designFile = paste0(root, "/design_file2.tsv")
   #outdir     = paste0(root, "/DOA_pairwise")
   #prefix     = "edger_otu"
   #pvalue     = 0.15
   #fdr        = 0.15
   #logfc      = 1.1
     
   print(paste0("[DEBUG] infile: ",infile))
   print(paste0("[DEBUG] outdir: ",outdir))
   print(paste0("[DEBUG] design_file: ",designFile))
 
   vColors = c( '#00FF00', '#FF8080', '#FF00FF', '#0000FF', '#808282', '#CCFFFF',
                '#CCFFCC', '#99CCFF', '#CC99FF', '#FFCC99', '#3366FF', '#33CCCC',
                '#99CC00', '#FF99CC', '#FFCC00')

   # load design file. First row of design file is labeled #group and contains comparisons that
   # are relevant as a group. 
   design_file_1 = read.table(designFile, sep="\t", comment.char="", nrows=1, colClasses="character", header=F, check.names=FALSE)
   factors = factor(as.character(design_file_1[1,2:ncol(design_file_1)]))
   myLevels = levels(factors)
   
   # Load design file and store header
   design_file = read.table(designFile, sep="\t", comment.char="", colClasses="character", header=T, check.names=FALSE)
     
   #load data and work the gene count table to make it compatible with edgeR
   tData = fread(infile, header=T, skip="#OTU ID", sep="\t", showProgress=FALSE)
   tData = data.frame(tData, check.names=FALSE)
   colnames(tData)[1] = "Symbol"
   tData2 = tData[, 2:ncol(tData)]
   rownames(tData2) = tData[ ,1]
   tDataTaxonomy = data.frame(tData2[,ncol(tData2)])
   rownames(tDataTaxonomy) = row.names(tData2)
   colnames(tDataTaxonomy)[1] = "taxonomy"
   tData2$taxonomy = NULL
   
   # We will store all DDA genes in this vector.
   DDA_list = vector()
   
   # First normalize all data together. Pairwise comparisons will be done later in inner loop
   #tmpGroup = c(rep("A", ncol(tData2)))
   #tData2 = tData2 + 1
   z <- DGEList(tData2, remove.zeros = TRUE)
   z <- calcNormFactors(z, method="TMM")
   w1 = cpm(z)
   # And generate normalized counts.
   #z2 = z
   #w = tData2
   #w1 = w
   #for(i in 1:ncol(w)){
   #  lib_size = z2$samples[c(names(w)[i]),]$lib.size
   #  norm_fac = z2$samples[c(names(w)[i]),]$norm.factors
   #  w1[,i] = w[,i] / (lib_size * norm_fac)
   #}
   
   #start analyses by groups (outer loop.)
   for(j in myLevels){
    
     DDA_list_curr_group = vector()
     samples_curr_group = vector()
     files = vector()
     files_normalized = vector()
     files_normalized_significant = vector()
          
     # For first level/group, prepare data frame
     curr_design_tmp = NULL
     curr_design_tmp = design_file[colnames(design_file) %in% j]
     curr_design = NULL
     curr_design$SampleID = design_file[,1]
     curr_design = data.frame(curr_design)
     curr_design = cbind(curr_design, curr_design_tmp)
     my_names = curr_design[1,]
     my_names$SampleID = NULL
     my_names = c("SampleID", as.character(my_names))
     curr_design = curr_design[-1,]
     colnames(curr_design) <- my_names
     
     for(i in 2:ncol(curr_design)){
             
       DDA_list_curr_design = vector()
       ncol(curr_design)
       experiment = colnames(curr_design)[i]
     
       currOutdir = paste0(outdir, "/design_group_", j, "/", experiment)
       dir.create(currOutdir, showWarnings=FALSE, recursive=TRUE)
             
       controls   = curr_design[curr_design[,i] == "1", 1]
       treatments = curr_design[curr_design[,i] == "2", 1]
       # Before contining, further parse controls and treatments to make sure they actually are in the tData2 (otu table) table.
       controls = controls[is.element(controls, colnames(tData2))]
       treatments = treatments[is.element(treatments, colnames(tData2))]
       
       sampleNames = c(controls,treatments)
       samples_curr_group = c(samples_curr_group, sampleNames)
       group1 = rep("Control", length(controls))
       group2 = rep("Treatment", length(treatments))
       group = c(group1, group2)
              
       # Extract columns corresponding to selected sample of current design.
       tData3 = tData2[,colnames(tData2) %in% sampleNames]
       #tData4 = tData3[sampleNames]
       tData4 = tData3[is.element(colnames(tData3), sampleNames)]
       tData4  = tData4[sampleNames]
     
       # classic edger
       #tData4 = tData4 + 1
       # reorder table according to design.
       y <- DGEList(tData4, group=group, remove.zeros=TRUE)
       ###### 
       keep <- rowSums(cpm(y)>1) >= 2
       y <- y[keep, , keep.lib.sizes=FALSE]
       #####
  
       ###  Continue Edger workflow.
       y1 <- calcNormFactors(y, method="TMM")
       y2 <- estimateCommonDisp(y1)
       y3 <- estimateTagwiseDisp(y2)
       et <- exactTest(y3)
       number_of_rows = nrow(et$table)

       # Work DEG to make it as a data frame to print as a tsv file.
       DEGs = topTags(et, n=number_of_rows)$table
       myColNames = c("symbol", colnames(DEGs))
       df <- NULL
       df$symbol = row.names(DEGs)
       df = data.frame(df)
       df$logFC = DEGs$logFC
       df$logCPM = DEGs$logCPM
       df$PValue = DEGs$PValue
       df$FDR = DEGs$FDR
       # Then filter for DEG having pvalue < 0.05
       df2 = df[df$PValue <= pvalue,]
       df2 = df2[df2$FDR <= fdr,]
       df2 = df2[abs(df2$logFC) >= logfc, ]
       df3 = merge(df2, tDataTaxonomy, by.x="symbol", by.y="row.names") 
       write.table(df3, file=paste0(currOutdir, "/", experiment, "_edger.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
       
       # Output normalized counts. i.e. Don't forget that lib.size and norm.factors come from variable z which are computed considering all samples.     
       x1 = tData4
       for(k in 1:ncol(tData4)){
         lib_size = y3$samples[c(names(tData4)[k]),]$lib.size
         norm_fac = y3$samples[c(names(tData4)[k]),]$norm.factors
         x1[,k] = tData4[,k] / (lib_size * norm_fac)
       } 
       write.table(x1, file=paste0(currOutdir, "/", experiment, "_normalized.tsv"), quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)
              
       # Then add file to list of vector
       files = c(files, paste0(currOutdir, "/", experiment, "_edger.tsv"))
       files_normalized = c(files_normalized, paste0(currOutdir, "/", experiment, "_normalized.tsv"))
          
       # Also write table containing normalized_significant
       df3 = x1[row.names(x1) %in% df2$symbol,]
       #write.table(df3, file=paste0(currOutdir, "/", experiment, "_normalized_significant.tsv"), quote=FALSE, sep="\t", row.names=TRUE, col.names=NA)
       #files_normalized_significant = c(files_normalized, paste0(currOutdir, "/", experiment, "_normalized_significant.tsv"))
       
       # Table for curr design only.
       samples_curr_design = unique(sampleNames)
       #DDA_table_curr_design = tData2[row.names(tData2) %in% df2$symbol, ] # raw counts (not normalized)
       DDA_table_curr_design = x1[row.names(x1) %in% df2$symbol, ] # raw counts (normalized)
       DDA_table_curr_design = DDA_table_curr_design[,samples_curr_design] # make sure to select only curr samples in the current design.
       
       DDA_df_curr_design = NULL
       DDA_df_curr_design = row.names(DDA_table_curr_design)
       DDA_df_curr_design = data.frame(DDA_df_curr_design)
       DDA_df_curr_design = cbind(DDA_df_curr_design, DDA_table_curr_design[,1:ncol(DDA_table_curr_design)])
       colnames(DDA_df_curr_design)[1] = "cluster"
       
       DDA_output_curr_design = paste0(outdir, "/design_group_", j, "/", experiment, "/DDA_normalized_significant.tsv")
       write.table(DDA_df_curr_design, file=(DDA_output_curr_design), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
      
       # Then write density outfiles.
       #if(nrow(df3) >= 2){
       #  density_plot_outfile_prefix = paste0(currOutdir, "/", experiment, "_normalized_significant_densityPlot")
       #  generate_density_plot(df3, controls, treatments, density_plot_outfile_prefix)
       #}
     }
     
     # Then once we have all files. Creat a table of all DEGs.
     # 1- Create a list of all genes diff. expressed.
     dfg = NULL
     genes = vector()
     genes_normalized = vector()
     for(x in files){
        print(x)
        table = read.table(x, sep="\t", header=T)
        genes = c(genes, table$symbol)
        DDA_list = c(DDA_list, table$symbol)
        DDA_list_curr_group = c(DDA_list_curr_group, table$symbol)
     }
     genes = unique(genes)
     dfg$symbol = genes
     dfg = data.frame(dfg)
   
     # 2- reloop through result table and merge results
     df_fc = dfg
     df_cpm = dfg
     df_pvalue = dfg
     df_fdr = dfg
     
     for(x in files){
       print(x) 
       compName = basename(file_path_sans_ext(x))
       table = read.table(x, sep="\t", header=T)
       
       #fc
       df_fc = merge(df_fc, table, by="symbol", all=TRUE)
       df_fc$logCPM = NULL
       df_fc$PValue = NULL
       df_fc$FDR = NULL
       names(df_fc)[names(df_fc) == 'logFC'] <- compName
       
       #cpm
       df_cpm = merge(df_cpm, table, by="symbol", all=TRUE)
       df_cpm$logFC = NULL
       df_cpm$PValue = NULL
       df_cpm$FDR = NULL
       names(df_cpm)[names(df_cpm) == 'logCPM'] <- compName
       
       #pvalue
       df_pvalue = merge(df_pvalue, table, by="symbol", all=TRUE)
       df_pvalue$logFC = NULL
       df_pvalue$logCPM = NULL
       df_pvalue$FDR = NULL
       names(df_pvalue)[names(df_pvalue) == 'PValue'] <- compName
       
       #fdr
       df_fdr = merge(df_fdr, table, by="symbol", all=TRUE)
       df_fdr$logFC = NULL
       df_fdr$logCPM = NULL
       df_fdr$PValue = NULL
       names(df_fdr)[names(df_fdr) == 'FDR'] <- compName     
     }
     
     curr_output_fc = paste0(outdir, "/design_group_", j, "/fc.tsv")
     curr_output_cpm = paste0(outdir, "/design_group_", j, "/cpm.tsv")
     curr_output_pvalue = paste0(outdir, "/design_group_", j, "/pvalue.tsv")
     curr_output_fdr = paste0(outdir, "/design_group_", j, "/fdr.tsv")
     write.table(df_fc, file=curr_output_fc, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
     write.table(df_cpm, file=curr_output_cpm, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
     write.table(df_pvalue, file=curr_output_pvalue, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
     write.table(df_fdr, file=curr_output_fdr, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
     
     # Table for curr group only.
     samples = unique(samples_curr_group)
     DDA_list_curr_group = unique(DDA_list_curr_group)
     #DDA_table_curr_group = tData2[row.names(tData2) %in% DDA_list_curr_group, ]
     DDA_table_curr_group = w1[row.names(w1) %in% DDA_list_curr_group, ] # raw counts (normalized)
     DDA_table_curr_group = DDA_table_curr_group[,samples] # make sure to select only curr samples in the current design.
     
     DDA_df_curr_group = NULL
     DDA_df_curr_group = row.names(DDA_table_curr_group)
     DDA_df_curr_group = data.frame(DDA_df_curr_group)
     DDA_df_curr_group = cbind(DDA_df_curr_group, DDA_table_curr_group[,1:ncol(DDA_table_curr_group)])
     colnames(DDA_df_curr_group)[1] = "cluster"
     
     DDA_output_curr_group = paste0(outdir, "/design_group_", j, "/DDA_normalized_significant.tsv")
     write.table(DDA_df_curr_group, file=DDA_output_curr_group, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
     
     
  }
  
  #Then create a normalized gene abundance table with only DDA genes in it. i.e. Collections of genes for all designs. 
  #Not sure it is useful as there are far too many genes in there.
  DDA_list = unique(DDA_list)
  DDA_table = w1[row.names(w1) %in% DDA_list, ]
  
  DDA_df = NULL
  DDA_df = row.names(DDA_table)
  DDA_df = data.frame(DDA_df)
  DDA_df = cbind(DDA_df, DDA_table[,1:ncol(DDA_table)])
  colnames(DDA_df)[1] = "cluster"
  
  DDA_output = paste0(outdir, "/DDA_normalized_significant.tsv")
  write.table(DDA_df, file=DDA_output, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
  
  print("[DEBUG] Done EdgeR R wrapper...")
} 

usage=function(errM) {
        cat("\nUsage : Rscript edger.R [option] <Value>\n")
        cat("       -i        : infile\n")
        cat("       -o        : outdir\n")
        cat("       -p        : pvalue cutoff\n")
        cat("       -l        : log fold change cutoff\n")
        cat("       -f        : fdr value cutoff\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 5) {
	usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
	if (ARG[i] == "-i") {
		infile=ARG[i+1]
	} else if (ARG[i] == "-o") {
		outdir=ARG[i+1]
	}else if (ARG[i] == "-d") {
		designFile=ARG[i+1]
	}else if (ARG[i] == "-p") {
	  pvalue=ARG[i+1]
   }else if (ARG[i] == "-f") {
	  fdr=ARG[i+1]
   }else if (ARG[i] == "-l") {
	  logfc=ARG[i+1]
   }
}

performEdgerPairwiseOTUs(infile, outdir, designFile, pvalue, fdr, logfc)

