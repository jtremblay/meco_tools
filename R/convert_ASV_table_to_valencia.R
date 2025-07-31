#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# Convert ASV table to ASV table for the VALENCIA software.
# The valencia output table is in .csv format.
# Author: Julien Tremblay - jtremblay514@gmail.com
convertASVTableToValencia <- function(infile, outfile) {
        
    library(data.table)
    library(dplyr)
    
    #outdir = "/project/6049199/projects/p110/FLO/"
    #setwd("/project/6049199/projects/p110/FLO/")
    #infile = "./amplicontagger/reads_12/feature_tables/feature_table_filtered_bacteriaArchaea_species.tsv"
    
    processTaxonomicLineages <- function(asv_table, summarize_lineage, tax_level="L6_only"){
      if(!is.logical(summarize_lineage)){
        stop("summarize_lineage variable needs to be a boolean (TRUE or FALSE)")
      }
      
      asv_table$TMP_ID = seq(1,nrow(asv_table),1)
      for(p in 1:nrow(asv_table)){
        y = strsplit(asv_table[[ "taxonomy" ]][p], ";", fixed=TRUE)[[1]]
        
        if(length(y) < 2){
          y[2] = "p__unknown"; y[3] = "c__unknown"; y[4] = "o__unknown"; y[5] = "f__unknown"; y[6] = "g__unknown"; y[7] = "s__unknown";
          
        }else if(length(y) < 3){
          y[3] = "c__unknown"; y[4] = "o__unknown"; y[5] = "f__unknown"; y[6] = "g__unknown"; y[7] = "s__unknown";
          
        }else if(length(y) < 4){
          y[4] = "o__unknown"; y[5] = "f__unknown"; y[6] = "g__unknown"; y[7] = "s__unknown";
          
        }else if(length(y) < 5){
          y[5] = "f__unknown"; y[6] = "g__unknown"; y[7] = "s__unknown";
          
        }else if(length(y) < 6){
          y[6] = "g__unknown"; y[7] = "s__unknown";
          
        }else if(length(y) < 7){
          y[7] = "s__unknown";
          
        }else if(length(y) < 8){
          
        }
        
        if(p == 1){
          tmp_df = data.frame(y)
          tmp_df = t(tmp_df)
          row.names(tmp_df) = asv_table[p ,c("TMP_ID")]
          
        }else{
          tmp_df = rbind(tmp_df, y)
          row.names(tmp_df)[p] = asv_table[p ,c("TMP_ID")]
        }
      }
      y2 = data.frame(tmp_df)
      
      colnames(y2) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      
      if(summarize_lineage == TRUE){
        if(tax_level == "L6_only"){
          tmp_taxons = data.frame(TMP_ID=row.names(y2), paste0(y2[,6]))
        }else if(tax_level == "L6"){
          tmp_taxons = data.frame(TMP_ID=row.names(y2), paste0(y2[,4], ";", y2[,6]))
        }else if(tax_level == "L5"){
          tmp_taxons = data.frame(TMP_ID=row.names(y2), paste0(y2[,2], ";", y2[,5]))
        }else if(tax_level == "L4"){
          tmp_taxons = data.frame(TMP_ID=row.names(y2), paste0(y2[,2], ";", y2[,4]))
        }else if(tax_level == "L3"){
          tmp_taxons = data.frame(TMP_ID=row.names(y2), paste0(y2[,2], ";", y2[,3]))
        }else if(tax_level == "L2"){
          tmp_taxons = data.frame(TMP_ID=row.names(y2), paste0(y2[,2]))
        }else if(tax_level == "L1"){
          tmp_taxons = y2[,1,drop=FALSE]
          tmp_taxons = data.frame(TMP_ID=row.names(y2), paste0(y2[,1]))
        }else if(tax_level == "L7"){
          tmp_taxons = data.frame(TMP_ID=row.names(y2), paste0(y2[,4], ";", y2[,6], ";", y2[,7]))
        }else if(tax_level == "L7_only"){
            tmp_taxons = data.frame(TMP_ID=row.names(y2), paste0(y2[,7]))
        }
          
        tmp_taxons = data.frame(tmp_taxons)
        colnames(tmp_taxons)[2] = "taxonomy"
        taxons = unique(tmp_taxons)
        tmp = merge(asv_table, tmp_taxons, by.x="TMP_ID", by.y="TMP_ID")
        tmp$taxonomy.x = NULL
        colnames(tmp)[ncol(tmp)] = "taxonomy"
        tmp$TMP_ID = NULL
        #tmp_colnames = colnames(tmp)[1:(ncol(tmp)-1)]
        #tmp_colnames = c("taxonomy", tmp_colnames)
        #tmp = tmp[,tmp_colnames]
        #asv_table = tmp
      }
      
      return(tmp)
    }
    
    getDeepestTaxonomicLineage <- function(asv_table){
      
      asv_table$TMP_ID = seq(1,nrow(asv_table),1)
      for(p in 1:nrow(asv_table)){
        y = strsplit(asv_table[[ "taxonomy" ]][p], ";", fixed=TRUE)[[1]]
        deepest = NULL
        if(length(y) == 1){
          deepest = y[1];
        }else if(length(y) == 2){
          deepest = y[2];
        }else if(length(y) == 3){
          deepest = y[3];
        }else if(length(y) == 4){
          deepest = y[4];
        }else if(length(y) == 5){
          deepest = y[5];
        }else if(length(y) == 6){
          deepest = y[6];
        }else if(length(y) == 7){
          deepest = y[7];
        }
        
        asv_table[p ,c("TMP_ID")] = deepest
      }
      
      asv_table$taxonomy = asv_table$TMP_ID; asv_table$TMP_ID = NULL;
      
      return(asv_table)
    }
    
    my_asv_table_file = infile
    my_asv_table = data.frame(fread(my_asv_table_file), check.names=F)
    head(my_asv_table)
    my_asv_table = getDeepestTaxonomicLineage(my_asv_table) 
    head(my_asv_table)
    my_asv_table$taxonomy = gsub("^s__", "", my_asv_table$taxonomy)
    # First split taxonomy
    colnames(my_asv_table)[1] = "feature_id"
    my_asv_table$feature_id = NULL
    #my_asv_table = my_asv_table[!my_asv_table$taxonomy %in% c("g__unknown"),]
    #tmp = my_asv_table[, c(1:10,381)]
    #head(tmp)
    #tmp$feature_id = NULL
    X2 = data.table::transpose(my_asv_table, keep.names="taxonomy", make.names="taxonomy")
    colnames(X2)[1] = "sampleID"
    X2$read_count = rowSums(X2[,c(2:ncol(X2))])
    X2 = X2[,c("sampleID", "read_count", colnames(X2)[!colnames(X2) %in% c("sampleID", "read_count")])]
    #write.table(X2, file=paste0(outdir, "/amplicontagger/reads_12/feature_tables/feature_table_for_valencia.csv"), row.names=F, sep=",", quote=F)
    write.table(X2, file=outfile, row.names=F, sep=",", quote=F)
}

usage=function(errM) {
    cat("\nUsage : Rscript convert_ASV_table_to_valencia.R [option] <Value>\n")
    cat("       -i        : ASV table in .tsv format\n")
    cat("       -o        : outfile (valencia tax table in .csv format\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 2) {
    usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-i") {
        infiles=ARG[i+1]
    } else if (ARG[i] == "-o") {
        outfile=ARG[i+1]
    }
}

usage=function(errM) {
  cat("\nUsage : Rscript mergeRawCounts.R [option] <Value>\n")
  cat("       -i        : List of files separated by a comma\n")
  cat("       -o        : outfile\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 2) {
  usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-i") {
    infile=ARG[i+1]
  } else if (ARG[i] == "-o") {
    outfile=ARG[i+1]
  }
}
convertASVTableToValencia(infile, outfile)