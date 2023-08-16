#!/usr/bin/env Rscript

# Function that takes input values and run multiple edgeR analyses based on blocks and treatments
# Julien Tremblay => julien.tremblay@nrc-cnrc.gc.ca
performEdgerGLM <- function(infile, outdir, mappingFile, pvalue, fdr, logfc, blockColumnLabels, treatmentColumnLabels) {
   library(edgeR)
   library(tools)
   library(data.table)

   options(stringsAsFactors = FALSE)
 
   ## For debugging 
   #root                  = "/home/jtrembla/Projects/test_R_block_effects/"
   #infile                = paste0(root, "merged_gene_abundance.tsv")
   #designFile            = paste0(root, "design_file_AcBarrie.tsv")
   #mappingFile           = paste0(root, "mapping_file.tsv")
   #outdir                = paste0(root, "/DDA_GLM/")
   #prefix                = "edger"
   #pvalue                = 0.05
   #fdr                   = 0.05
   #logfc                 = 1.5
   #blockColumnLabels     = "Block"
   #treatmentColumnLabels = "Variety"

   logfc = as.numeric(logfc)
   fdr = as.numeric(fdr)
   pvalue = as.numeric(pvalue)
   
   dir.create(outdir, showWarnings = FALSE)
    
   print(paste0("[DEBUG] infile: ", infile))
   print(paste0("[DEBUG] outdir: ", outdir))
   print(paste0("[DEBUG] mapping file: ", mappingFile))
 
   vColors = c( '#00FF00', '#FF8080', '#FF00FF', '#0000FF', '#808282', '#CCFFFF',
                '#CCFFCC', '#99CCFF', '#CC99FF', '#FFCC99', '#3366FF', '#33CCCC',
                '#99CC00', '#FF99CC', '#FFCC00')
 
   # Load mapping file. Especially to find blocks
   mapping_file = read.table(mappingFile, sep="\t", colClasses = "character")
   header = read.table(mappingFile, sep="\t", nrows=1, comment.char="")
   colnames(mapping_file) = header
   row.names(mapping_file) = mapping_file[,1]
   mapping_file[,1] = NULL
   mapping_file_original_names = colnames(mapping_file)
 
   # Check if supplied blockColumnLabels and treatmentColumns labels are present in mapping file.
   # Split blockColumnLabels and treatmentColumnsLabels by ','
   block_columns <- strsplit(blockColumnLabels, ",")[[1]]
   treatment_columns <- strsplit(treatmentColumnLabels, ",")[[1]]
   for(test_block in block_columns){
      print(paste0("test_block:", test_block))
      if(test_block != "no"){
         if(!(test_block %in% colnames(mapping_file))){
            stop(paste0("One of your Block -b <", test_block, "> variable is not in the abundance matrix...\n"))
         }
      }
   }
   for(test_treatment in treatment_columns){
      if(!(test_treatment %in% colnames(mapping_file))){
         stop(paste0("One of your Treatment -t <", test_treatment, "> variable is not in the abundance matrix...\n"))
      }
   }
    

   #load data and work the gene count table to make it compatible with edgeR
   tData = fread(infile, header=T, skip=0, sep="\t", showProgress=FALSE)
   tData = data.frame(tData, check.names=FALSE)
   colnames(tData)[1] = "Symbol"
   tData2 = tData[, 2:ncol(tData)]
   rownames(tData2) <- tData[ ,1]
   
   # Only keep rows that are actually in the mapping file.
   # So here order of sample IDs in mapping file is the same as sample IDs in tData2 (count table)
   #tData2 = tData2[row.names(mapping_file)]
   tData2 = tData2[is.element(colnames(tData2), row.names(mapping_file))]
   
   # Also, in mapping_file : only keep samples that are present in otu table.
   if(ncol(mapping_file) == 1){ # fix when mapping file is a df containing only one column.
      mapping_file$dummy = "BLABLA"
      mapping_file = mapping_file[is.element(row.names(mapping_file), colnames(tData2)),]
      mapping_file$dummy = NULL
   
   }else{
      mapping_file = mapping_file[is.element(row.names(mapping_file), colnames(tData2)),]
   }
   print(head(mapping_file))
   #print("row.names(mapping_file)")
   #print(row.names(mapping_file))
   #print("colnames(tData2)")
   #print(colnames(tData2))
   

   #print(paste0("block_columns:", block_columns))
   #print(paste0("treatment_columns:", treatment_columns))
   tData2 = tData2 + 1
   # Their could be multiple 'treatment' columns.
   for(curr_treatment_column in treatment_columns){
      print(paste0("Computing GLM for : ", curr_treatment_column))
      groups = mapping_file[[curr_treatment_column]]
      y <- DGEList(counts=tData2, remove.zeros=TRUE)

      # We filter out lowly expressed genes using the following commands.
      # EdgeR recommandation.
      keep <- rowSums(cpm(y)>1) >= 2
      y <- y[keep, , keep.lib.sizes=FALSE]
      
      #treatment_list = sort(mapping_file[[curr_treatment_column]]) 
      #treatments = factor(treatment_list)
      #unique_treatments = unique(treatments)
      #unique_treatments = as.character(unique_treatments)

      treatment_df = data.frame(row.names(mapping_file), mapping_file[[curr_treatment_column]])
      colnames(treatment_df) = c("sampleID", "curr_treatment_column")
      treatment_df = treatment_df[match(colnames(y$counts), treatment_df$sampleID), ]
      treatment_list = as.character(treatment_df$curr_treatment_column)
      treatments = factor(treatment_list)
      unique_treatments = levels(treatments)
      #print("treatments:")
      #print(treatments)
     
      #print(paste0("ncol(tData2): ", ncol(tData2)))
      #print(paste0("nrow(mapping_file): ", nrow(mapping_file)))
 
      # Then extract appropriate variables
      # Then extract appropriate variables
      if(identical(block_columns, character(0)) | block_columns == "no" ){ # means no blocks.

         outdir = paste0(outdir, "/contrasts")
         dir.create(outdir, showWarnings = FALSE)

         design <- model.matrix(~0+treatments, data=y$samples)
         
         y <- calcNormFactors(y, method="TMM")
         y <- estimateGLMCommonDisp(y,design)
         y <- estimateGLMTrendedDisp(y,design)
         y <- estimateGLMTagwiseDisp(y,design)
         fit <- glmFit(y,design)
         
         # Then for each treatment, Let's loop through treatment_list by releveling with each treatments as they come in...         
         for(i in 1:(length(unique_treatments)-1)){
            curr_control = unique_treatments[i]
            print(paste0("Current control: ", curr_control))
            
            # Here we'll proceed by contrast. for instance if we have 7 columns (i.e. treatments)
            #contrast=c(-1,1,0,0,0,0,0) # -1=control, 1=treatment, and 0=leave out.. Will compare 2 column vs 1 column
            
            #unique_treatments2 = unique_treatments[(i+1):(length(unique_treatments)-0)]
            #unique_treatments2 = unique_treatments
            
            #Then from our unique_treatments id, find indices of these columns in the design file.
            
            # Make directory for current ref/control
            curr_outdir = paste0(outdir, '/')#, as.character(treatments[1]))
            dir.create(curr_outdir, showWarnings=FALSE, recursive=TRUE)
            for(j in (i+1):length(unique_treatments)){

               curr_contrast = rep(0, length(unique_treatments))
               curr_contrast[i] = -1
               curr_contrast[j] = 1
               curr_treatment = unique_treatments[j]
               
               lrt <- glmLRT(fit, contrast=curr_contrast)
               DEGs = topTags(lrt, n=nrow(tData2))$table
               curr_col_names = c("symbol", colnames(DEGs))
               DEGs1 <- NULL
               DEGs1$symbol = row.names(DEGs)
               DEGs1 = data.frame(DEGs1)
               DEGs1$logFC = DEGs$logFC
               DEGs1$logCPM = DEGs$logCPM
               DEGs1$PValue = DEGs$PValue
               DEGs1$FDR = DEGs$FDR
               
               # Then filter for DEG having pvalue < 0.05
               DEGs1 = DEGs1[DEGs1$PValue <= pvalue,]
               DEGs1 = DEGs1[DEGs1$FDR <= fdr,]
               DEGs1 = DEGs1[abs(DEGs1$logFC) >= logfc, ]
               DEGs1 = DEGs1[complete.cases(DEGs1),]
               print(paste0(curr_outdir, "/", curr_treatment,"_vs_", curr_control, "_edger.tsv"))
               
               write.table(DEGs1, file=paste0(curr_outdir, "/", curr_treatment,"_vs_", curr_control, "_edger.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
               
               # Then do this for invert contrast (i.e. to generate all possibilities.)               
               lrt <- glmLRT(fit, contrast=(curr_contrast * -1))
               DEGs = topTags(lrt, n=nrow(tData2))$table
               curr_col_names = c("symbol", colnames(DEGs))
               DEGs1 <- NULL
               DEGs1$symbol = row.names(DEGs)
               DEGs1 = data.frame(DEGs1)
               DEGs1$logFC = DEGs$logFC
               DEGs1$logCPM = DEGs$logCPM
               DEGs1$PValue = DEGs$PValue
               DEGs1$FDR = DEGs$FDR
               
               # Then filter for DEG having pvalue < 0.05
               DEGs1 = DEGs1[DEGs1$PValue <= pvalue,]
               DEGs1 = DEGs1[DEGs1$FDR <= fdr,]
               DEGs1 = DEGs1[abs(DEGs1$logFC) >= logfc, ]
               DEGs1 = DEGs1[complete.cases(DEGs1),]
               print(paste0(curr_outdir, "/", curr_control,"_vs_", curr_treatment, "_edger.tsv"))
               
               # Assign taxonomy to each DOA OTUs.
               write.table(DEGs1, file=paste0(curr_outdir, "/", curr_control,"_vs_", curr_treatment, "_edger.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)            
            }
         }
      }else{
         outdir = paste0(outdir, "/intercept")
         dir.create(outdir, showWarnings = FALSE)
         
         treatment_df2 = data.frame(row.names(mapping_file), mapping_file[[curr_treatment_column]], mapping_file[[block_columns]])
         colnames(treatment_df2) = c("sampleID", "curr_treatment_column", "block")
         treatment_df2 = treatment_df2[match(colnames(y$counts), treatment_df2$sampleID), ]
         block_list = as.character(treatment_df2$block)
         block = factor(block_list)
         unique_blocks = levels(block)
         
         design <- model.matrix(~block+treatments, data=y$samples)   
      
         y <- calcNormFactors(y, method="TMM")
         y <- estimateGLMCommonDisp(y,design)
         y <- estimateGLMTrendedDisp(y,design)
         y <- estimateGLMTagwiseDisp(y,design)
         fit <- glmFit(y,design)

         # Then for each treatment, Let's loop through treatment_list by releveling with each treatments as they come in...
         for(i in 1:(length(unique_treatments)-1)){
            curr_control = unique_treatments[i]
            print(paste0("Current control: ", curr_control))
               
            # Here we have n blocks, so compute DEGs by adjusting for block differences; so basically by leaving
            # out columns 1 to n which corresponds to the n blocks.
            # So we first need to know what variables needs to be considered as blocks. \
            
            unique_treatments2 = unique_treatments[(i+1):(length(unique_treatments)-0)]
            #Then from our unique_treatments id, find indices of these columns in the design file.
            
            # Make directory for current ref/control
            curr_outdir = paste0(outdir, '/')#, as.character(treatments[1]))
            dir.create(curr_outdir, showWarnings=FALSE, recursive=TRUE)
            for(j in 1:length(unique_treatments2)){
               full_design_name = paste0('treatments', unique_treatments2[j])
               index = which(colnames(design)==full_design_name)
               lrt <- glmLRT(fit, coef=index)
               DEGs = topTags(lrt, n=nrow(tData2))$table
               curr_col_names = c("symbol", colnames(DEGs))
               DEGs1 <- NULL
               DEGs1$symbol = row.names(DEGs)
               DEGs1 = data.frame(DEGs1)
               DEGs1$logFC = DEGs$logFC
               DEGs1$logCPM = DEGs$logCPM
               DEGs1$PValue = DEGs$PValue
               DEGs1$FDR = DEGs$FDR
               
               # Then filter for DEG having pvalue < 0.05
               DEGs1 = DEGs1[DEGs1$PValue <= pvalue,]
               DEGs1 = DEGs1[DEGs1$FDR <= fdr,]
               DEGs1 = DEGs1[abs(DEGs1$logFC) >= logfc, ]
               DEGs1 = DEGs1[complete.cases(DEGs1),]
               print(paste0(curr_outdir, "/", unique_treatments2[j],"_vs_", curr_control, "_edger.tsv"))
               
               # Assign taxonomy to each DOA OTUs.
               write.table(DEGs1, file=paste0(curr_outdir, "/", unique_treatments2[j],"_vs_", curr_control, "_edger.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)            
            }
         }
      }
   }
   print("[DEBUG] Done EdgeR R wrapper...")
} 

usage=function(errM) {
        cat("\nUsage : Rscript edger.R [option] <Value>\n")
        cat("       -i        : infile\n")
        cat("       -o        : outdir\n")
        cat("       -m        : mapping file\n")
        cat("       -p        : pvalue cutoff\n")
        cat("       -l        : log fold change cutoff\n")
        cat("       -f        : fdr value cutoff\n")
        cat("       -b        : Block column labels (seperated by a ',' if many.\n")
        cat("       -t        : Treatment column labels (seperated by a ',' if many.\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 7) {
	usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
	if(ARG[i] == "-i") {
		infile=ARG[i+1]
	}else if (ARG[i] == "-o") {
		outdir=ARG[i+1]
	}else if (ARG[i] == "-m") {
		mapping_file=ARG[i+1]
	}else if (ARG[i] == "-p") {
	  pvalue=ARG[i+1]
   }else if (ARG[i] == "-f") {
	  fdr=ARG[i+1]
   }else if (ARG[i] == "-l") {
	  logfc=ARG[i+1]
   }else if (ARG[i] == "-b") {
      blocks=ARG[i+1]
   }else if (ARG[i] == "-t") {
      treatments=ARG[i+1]
   }
}

performEdgerGLM(infile, outdir, mapping_file, pvalue, fdr, logfc, blocks, treatments)

