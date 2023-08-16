#!/usr/bin/env Rscript

# Function that takes a infile: (Qiime's tax summary spreadsheets) and an outdir where to
# To write plot files. Julien Tremblay => julien.tremblay@nrc-cnrc.gc.ca
performEdger <- function(infile, outdir, designFile) {
   library(edgeR)
   library(tools)

   options(stringsAsFactors = FALSE)
   
   #infile = "~/Projects/metagenome_stats/merged_gene_abundance_50000.tsv"
   #outdir = "~/Projects/metagenome_stats"
   prefix = "edger"
   #designFile = "~/Projects/metagenome_stats/design_file_parsed.tsv"
  
   print(paste0("[DEBUG] infile: ",infile))
   print(paste0("[DEBUG] outdir: ",outdir))
   print(paste0("[DEBUG] design_file: ",designFile))
 
   vColors = c( '#00FF00', '#FF8080', '#FF00FF', '#0000FF', '#808282', '#CCFFFF',
                '#CCFFCC', '#99CCFF', '#CC99FF', '#FFCC99', '#3366FF', '#33CCCC',
                '#99CC00', '#FF99CC', '#FFCC00')

   #load design file. First row of design file is labeled #group and contains comparisons that
   #are relevant as a group. 
   design_file_1 = read.table(designFile, sep="\t", comment.char="", nrows=1, colClasses="character", header=F, check.names=FALSE)
   factors = factor(as.numeric(design_file_1[1,]))
   myLevels = levels(factors)
   
   # Load design file and store header
   design_file = read.table(designFile, sep="\t", comment.char="", colClasses="character", header=T, check.names=FALSE)
   #header = read.table(designFile, sep="\t", nrows=1, comment.char="", header=FALSE)
   #colnames(design_file) = header
   
   #load data and work the gene count table to make it compatible with edgeR
   tData <- read.csv(file=infile, header=T, sep="\t", check.names=FALSE)
   colnames(tData)[1] = "Symbol"
   tData2 = tData[, 2:ncol(tData)]
   rownames(tData2) <- tData[ ,1]
   
   #start analyses by groups (outer loop.)
   for(j in myLevels){
    
     files = vector()
          
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
             
       ncol(curr_design)
       experiment = colnames(curr_design)[i]
       #experiments = strsplit(experiment, split = "_vs_", perl = TRUE)
       #control_name = experiments[[1]][1]
       #treatment_name = experiments[[1]][2]
     
       currOutdir = paste0(outdir, "/design_group_", j, "/", experiment)
       dir.create(currOutdir, showWarnings=FALSE, recursive=TRUE)
       #setwd(currOutdir)
       
       controls = curr_design[design_file[,i] == "1", 1]
       treatments = curr_design[design_file[,i] == "2", 1]
       sampleNames = c(controls,treatments)
       group1 = rep("Control", length(controls))
       group2 = rep("Treatment", length(treatments))
       group = c(group1, group2)
     
       # Extract columns corresponding to selected sample of current design.
       tData3 = tData2[,colnames(tData2) %in% sampleNames]
       tData4 = tData3[sampleNames]
        
       # classic edger
       #tData5 = tData4[colnames(tData2) %in% mapping_file_parsed[,1]]
       y <- DGEList(tData4, group=group)
       y <- calcNormFactors(y)
       y <- estimateCommonDisp(y)
       y <- estimateTagwiseDisp(y)
       et <- exactTest(y)
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
       df2 = df[df$PValue < 0.01,]
       write.table(df2, file=paste0(currOutdir, "/", experiment, ".tsv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
     
       # Then add file to list of vector
       files = c(files, paste0(currOutdir, "/", experiment, ".tsv"))
          
       # glm edger
       #design <- model.matrix(~group)
       #z <- estimateGLMCommonDisp(y,design)
       #z <- estimateGLMTrendedDisp(y,design)
       #z <- estimateGLMTagwiseDisp(y,design)
       #fit <- glmFit(y,design)
       #lrt <- glmLRT(fit,coef=2)
       #topTags(lrt)
     }
     
     # Then once we have all files. Creat a table of all DEGs.
     # 1- Create a list of all genes diff. expressed.
     dfg = NULL
     genes = vector()
     for(x in files){
        print(x) 
        table = read.table(x, sep="\t", header=T)
        genes = c(genes, table$symbol)
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
       df_fdr = merge(df_pvalue, table, by="symbol", all=TRUE)
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
  }
  
  print("[DEBUG] Done EdgeR R wrapper...")
} 

usage=function(errM) {
        cat("\nUsage : Rscript taxBarplotWithMappingFile.R [option] <Value>\n")
        cat("       -i        : infile\n")
        cat("       -o        : outdir\n")
        cat("       -p        : prefix\n")
        cat("       -m        : design_file\n")
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
	}else if (ARG[i] == "-d") {
		designFile=ARG[i+1]
	}#else if (ARG[i] == "-p") {
	 # prefix=ARG[i+1]
	#}
}

performEdger(infile, outdir, designFile)

