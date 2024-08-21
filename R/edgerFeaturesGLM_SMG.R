#!/usr/bin/env Rscript

# Function that takes input values and run multiple edgeR analyses based on blocks and treatments
# National Research Council Canada
# Genomics and Microbiomes
# Julien Tremblay => julien.tremblay@nrc-cnrc.gc.ca
performEdgerGLMFeatures <- function(infile, outdir, mappingFile, pvalue, fdr, logfc, blockColumnLabels, treatmentColumnLabels, cutoff=150) {
   library(edgeR)
   library(tools)
   library(data.table)
   library(statmod)

   options(stringsAsFactors = FALSE)  
   removeSumZero = function(M) M[, colSums(abs(M)) != 0] 
 
   ## For debugging 
   #root                  = "/project/6004719/projects/GROW/Mesocosm_exp1/export/"
   #infile                = paste0(root, "gene_abundance/merged_gene_abundance.tsv")
   #mappingFile           = paste0(root, "mapping_file_DEG_CM1.tsv")
   #outdir                = paste0(root, "/DDA_GLM/")
   #pvalue                = 0.05
   #fdr                   = 0.1
   #logfc                 = 1.0
   #blockColumnLabels     = "no"
   #treatmentColumnLabels = "Treatment1,Treatment2"
   #cutoff = 150
   cutoff = as.numeric(cutoff)
   pvalue = as.numeric(pvalue)
   fdr = as.numeric(fdr)
   logfc = as.numeric(logfc)
      
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
  
   print(head(mapping_file))

   # Check if supplied blockColumnLabels and treatmentColumns labels are present in mapping file.
   # Split blockColumnLabels and treatmentColumnsLabels by ','
   block_columns <- strsplit(blockColumnLabels, ",")[[1]]
   treatment_columns <- strsplit(treatmentColumnLabels, ",")[[1]]
   
   if(identical(block_columns, character(0))){
      print("No block columns found...")
   }else{
      for(test_block in block_columns){
         if(test_block != "no"){
            if(!(test_block %in% colnames(mapping_file))){
               stop(paste0("One of your Block -b <", test_block, "> variable is not in the mapping file...\n"))
            }
         }
      }
   }
   for(test_treatment in treatment_columns){
      print(paste0("test_treatment:", test_treatment))
      if(!(test_treatment %in% colnames(mapping_file))){
         stop(paste0("One of your Treatment -t <", test_treatment, "> variable is not in the mapping file...\n"))
      }
   }

   #26423794      191
   #load data and work the gene count table to make it compatible with edgeR
   my_samples = unique(row.names(mapping_file))
   header = as.character(read.table(infile, header=F, nrows=1))
   tData = fread(infile, header=T, sep="\t", showProgress=T, select=c(header[1], my_samples))
   nrow(tData)
   ncol(tData)
   tData = data.frame(tData, check.names=FALSE)
   colnames(tData)[1] = "Symbol"
   row.names(tData) = tData$Symbol; tData$Symbol = NULL;
   
   # Only keep rows that are actually in the mapping file.
   # So here order of sample IDs in mapping file is the same as sample IDs in tData2 (count table)
   tData = tData[rowSums(tData) >= cutoff,]
   tData = tData[is.element(colnames(tData), row.names(mapping_file))]
   print(head(tData))
   print(dim(tData))
   # Also, in mapping_file : only keep samples that are present in feature table.
   mapping_file = mapping_file[is.element(row.names(mapping_file), colnames(tData)),]
   
    # Their could be multiple 'treatment' columns.
   tData = removeSumZero(tData) # Remove samples having zero counts!
   #tData2$filt<-apply(tData2, 1, function(x) sum(x>3)) #you can change the greater than to less than if you want to invert the count.
   #tData2 = tData2[tData2$filt>=3,] #the 2 is made up by me for the case of wanting 2 or more columns that are .005 or greater.  Change the 2 for your needs
   #tData2$filt<-NULL
   tData = tData + 1
   
   for(curr_treatment_column in treatment_columns){
      #curr_treatment_column = treatment_columns[2]
      print(paste0("Computing GLM for : ", curr_treatment_column))
      groups = mapping_file[[curr_treatment_column]]
      print("y before:") 
      y <- DGEList(counts=tData, remove.zeros=TRUE, group=groups)
      print("y after:") 
      print(y)
     
      # We filter out lowly expressed genes using the following commands.
      # EdgeR recommandation. But apparently, don 't do this for Features DA assessment?
      ################
      #keep <- rowSums(cpm(y)>1) >= 2
      #y <- y[keep, , keep.lib.sizes=FALSE]
      ################
      
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
      if(identical(block_columns, character(0)) | block_columns == "no" ){ # means no blocks.

         outdir1 = paste0(outdir, "/contrasts")
         dir.create(outdir, showWarnings = FALSE)
         
         design <- model.matrix(~0+treatments, data=y$samples)
         y <- calcNormFactors(y, method="RLE")
         # GLM pipeline
         #y <- estimateDisp(y,design)
         y <- estimateGLMCommonDisp(y,design)
         #y <- estimateGLMTrendedDisp(y,design)#, df=2)
         y <- estimateGLMTagwiseDisp(y,design)#, prior.df=10)
         # Normal pipeline
         #y <- calcNormFactors(y, method="RLE")
         #y <- estimateCommonDisp(y,design)
         #y <- estimateTagwiseDisp(y,design, prior.df=10)
         
         #y <- estimateDisp(y,design,robust=TRUE)#, prior.df=10)
         fit = glmFit(y,design)
         #fit = glmFit(y,design, prior.count=1)
         #fit <- glmQLFit(y,design, robust=TRUE, prior.count=2)#, winsor.tail.p=c(0.05, 0.1))
         
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
            curr_outdir = paste0(outdir1, '/')#, as.character(treatments[1]))
            dir.create(curr_outdir, showWarnings=FALSE, recursive=TRUE)
            for(j in (i+1):length(unique_treatments)){

               curr_contrast = rep(0, length(unique_treatments))
               curr_contrast[i] = -1
               curr_contrast[j] = 1
               curr_treatment = unique_treatments[j]
               
               lrt <- glmLRT(fit, contrast=curr_contrast)
               #lrt <- glmQLFTest(fit, contrast=curr_contrast)
               DEGs = topTags(lrt, n=nrow(tData))$table
      
               curr_col_names = c("symbol", colnames(DEGs))
               DEGs1 <- NULL
               DEGs1$symbol = row.names(DEGs)
               DEGs1 = data.frame(DEGs1)
               DEGs1$logFC = DEGs$logFC
               DEGs1$logCPM = DEGs$logCPM
               DEGs1$PValue = DEGs$PValue
               DEGs1$FDR = DEGs$FDR
               
               # Then filter for DEG having pvalue < 0.05
               DEGs1 = DEGs1[DEGs1$PValue <= pvalue, ]
               DEGs1 = DEGs1[DEGs1$FDR <= fdr, ]
               DEGs1 = DEGs1[abs(DEGs1$logFC) >= logfc, ]
               #DEGs1 = DEGs1[complete.cases(DEGs1), ]
               print(paste0(curr_outdir, "/", curr_treatment,"_vs_", curr_control, "_edger.tsv"))
               
               # Assign taxonomy to each DOA Features.
               #DEGs2 = merge(DEGs1, tDataTaxonomy, by.x="symbol", by.y="row.names")            
               write.table(DEGs1, file=paste0(curr_outdir, "/", curr_treatment,"_vs_", curr_control, "_edger.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
               
               # Then do this for invert contrast (i.e. to generate all possibilities.)               
               lrt <- glmLRT(fit, contrast=(curr_contrast * -1))
               #lrt <- glmQLFTest(fit, contrast=(curr_contrast * -1))
               DEGs = topTags(lrt, n=nrow(tData))$table
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
               #DEGs1 = DEGs1[complete.cases(DEGs1),]
               print(paste0(curr_outdir, "/", curr_control,"_vs_", curr_treatment, "_edger.tsv"))
               
               # Assign taxonomy to each DOA Features.
               #DEGs2 = merge(DEGs1, tDataTaxonomy, by.x="symbol", by.y="row.names")            
               write.table(DEGs1, file=paste0(curr_outdir, "/", curr_control,"_vs_", curr_treatment, "_edger.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)            
            }
         }
   
      # Else, do block design with intercept kind of design.
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
               lrt <- glmLRT(fit, coef=index, robuts=)
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
               #DEGs1 = DEGs1[complete.cases(DEGs1),]
               print(paste0(curr_outdir, "/", unique_treatments2[j],"_vs_", curr_control, "_edger.tsv"))
               
               # Assign taxonomy to each DOA Features.
               DEGs2 = merge(DEGs1, tDataTaxonomy, by.x="symbol", by.y="row.names")            
               write.table(DEGs2, file=paste0(curr_outdir, "/", unique_treatments2[j],"_vs_", curr_control, "_edger.tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)            
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

performEdgerGLMFeatures(infile, outdir, mapping_file, pvalue, fdr, logfc, blocks, treatments)

