#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# Generate OTU table from shotgun metagenomics
# reads abundance data.
# National Research Council - Biomonitoring
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
generateOTUTable <- function(gene_abun_filename, taxonomy_filename, outfile) {
   taxonomy = read.table(taxonomy_filename, sep="\t", header=T, comment.char="")
   colnames(taxonomy)[1] = "gene_id"

   #print(head(taxonomy))

   # exit if taxonomy table is empty.
   if(nrow(taxonomy) == 0){
      final_df = NULL
      write.table(final_df, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

   }else{
      # only if taxonomy is empty, read abundance file.
      gene_abun = read.table(gene_abun_filename, sep="\t", header=T, skip=0, comment.char="", check.names=FALSE)
      colnames(gene_abun)[1] = "gene_id"

      #print(head(gene_abun))

      ## Parse taxonomy
      taxonomy$kingdom = paste0("k__", taxonomy$kingdom)
      taxonomy$phylum = paste0("p__", taxonomy$phylum)
      taxonomy$class = paste0("c__", taxonomy$class)
      taxonomy$order = paste0("o__", taxonomy$order)
      taxonomy$family = paste0("f__", taxonomy$family)
      taxonomy$genus = paste0("g__", taxonomy$genus)
      taxonomy$species = paste0("s__", taxonomy$species)
      taxonomy$lineage = paste0(taxonomy$kingdom, ";", taxonomy$phylum, ";", taxonomy$class, ";", taxonomy$order, ";", taxonomy$family, ";", taxonomy$genus, ";", taxonomy$species)
      
      ## Create new table and add gene abundance to it.
      df = NULL
      df$gene_id = taxonomy$gene_id
      #g = head(gene_abun, 10000)

      df_merged = merge(df, gene_abun, by.x="gene_id", by.y="gene_id", all=FALSE)
      gene_id = df_merged$gene_id
      df_merged = round(df_merged[,-1,drop=FALSE], digits=0)
      df_merged = cbind(gene_id, df_merged)
      final_df = merge(df_merged, taxonomy[, c("gene_id","lineage")], by.x="gene_id", by.y="gene_id", all=FALSE)
      names(final_df)[1] = "#OTU ID"
      names(final_df)[ncol(final_df)] = "taxonomy"
      #print(head(final_df))
      # Convert chr to numeric:
      for(i in 2:(ncol(final_df)-1)){
         final_df[,i] = as.numeric(final_df[,i])
      }
      final_df = final_df[rowSums(final_df[,2:(ncol(final_df)-1),drop=FALSE]) != 0, , drop=FALSE] 
   
      write.table(final_df, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
   }
}

usage=function(errM) {
  cat("\nUsage : RscriptgenerateOTUTable.R .R [option] <Value>\n")
  cat("       -a        : gene abundance file\n")
  cat("       -t        : taxonomy table\n")
  cat("       -o        : outfile\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 3) {
  usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-a") {
    gene_abun=ARG[i+1]
  } else if (ARG[i] == "-t") {
    taxonomy=ARG[i+1]
  }else if (ARG[i] == "-o") {
    outfile=ARG[i+1]
  }
}

generateOTUTable(gene_abun, taxonomy, outfile)
