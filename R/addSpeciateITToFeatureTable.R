#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

# Add speciateIT to feature table.
# Author: Julien Tremblay - jtremblay514@gmail.com
addSpeciateITToFatureTable = function(infile, speciateIT_file, cutoff_pprob=0.8, outfile) {

    library(data.table)
    
    #infile = "./amplicontagger/reads_12/feature_tables/feature_table_filtered_bacteriaArchaea.tsv"
    #speciateIT_file = "./amplicontagger/reads_12/speciateIT/MC_order7_results.txt"
    #outfile = "./amplicontagger/reads_12/feature_tables/feature_table_filtered_bacteriaArchaea_species.tsv"
    #cutoff_pprob = 0.8

    df = data.frame(fread(infile), check.names=F)
    head(df)
      
    speciateIT = data.frame(fread(speciateIT_file), check.names=F)
    colnames(speciateIT) = c("FEATURE_ID", "species", "post_prob", "number_of_decisions")
    head(speciateIT) 
    
    df2 = merge(df, speciateIT, by.x="#FEATURE_ID", by.y="FEATURE_ID", all.x=T, all.y=F)
    head(df2)
    idx = which(df2$post_prob >= cutoff_pprob & grepl(";g__", df2$taxonomy))
    df2$taxonomy2 = df$taxonomy
    df2$taxonomy2[idx] = paste0(df2$taxonomy[idx], ";s__", df2$species[idx])
    df2$taxonomy2 = gsub(";;", ";", df2$taxonomy2)
    df2$taxonomy = NULL
    df2$species = NULL
    df2$post_prob = NULL    
    df2$number_of_decisions = NULL
    colnames(df2)[(ncol(df2))] = "taxonomy"
    
    fwrite(df2, outfile, sep="\t", quote=F, row.names=F)
}

usage=function(errM) {
    cat("\nUsage : Rscript addSpeciateITToFatureTable.R [option] <Value>\n")
    cat("       -i        : feature table (classified to the genus level (L6)\n")
    cat("       -c        : cutoff_pprob. Default = 0.8\n")
    cat("       -s        : speciateIT results file\n")
    cat("       -o        : outfile\n")
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
        outfile=ARG[i+1]
    } else if (ARG[i] == "-s") {
        speciateIT_file=ARG[i+1]
    } else if (ARG[i] == "-c") {
        cutoff_pprob=ARG[i+1]
    }
}

addSpeciateITToFatureTable(infile, speciateIT_file, cutoff_pprob, outfile)