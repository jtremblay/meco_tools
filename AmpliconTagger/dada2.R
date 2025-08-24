#!/usr/bin/env Rscript

# R wrapper for dada2.
# National Research Council Canada
# Genomics and Microbiomes
# Author: Julien Tremblay - julien.tremblay@nrc-cnrc.gc.ca
# type can be : spe (short paired-end), sse (short single-end) or lse (long single-end).
doDada2 <- function(indir, outdir, type, num_threads,
                    trimLeft.a, maxEE.a, truncQ.a, minQ.a, min_overlap) {

   library(dada2);packageVersion("dada2")
   #library(Biostrings)

   # Validation
   if(maxEE.a == -1){ maxEE.a = Inf }else{ maxEE.a = as.numeric(maxEE.a) }
   trimLeft.a  = as.numeric(trimLeft.a)
   truncQ.a    = as.numeric(truncQ.a)
   minQ.a      = as.numeric(minQ.a)
   min_overlap = as.numeric(min_overlap)
   
   print(paste0("Running DADA2 with the following arguments:",
             "trimLeft=", trimLeft.a," ",
             "maxEE=", maxEE.a," ",
             "truncQ=", truncQ.a," ",
             "num_threads=", num_threads," ",
             "type=", type," ",
             "indir=", indir," ",
             "min_overlap=", min_overlap, " ",
             "outdir=", outdir," "
            )
      
   )
  
   num_threads = as.integer(num_threads)

   dir.create(outdir)
   fns <- list.files(indir, pattern="fastq", full.names=TRUE)
   print("Found fastq files:")
   print(fns)
   rc <- dada2:::rc
   
   if(type == "lse"){
      
      lens.fn <- lapply(fns, function(fn) nchar(getSequences(fn)))
      lens <- do.call(c, lens.fn)
     
      png(filename=paste0(outdir, "/histogram.png"))
      hist(lens, 100, xlim = range(500,1600))
      dev.off()
      pdf(file=paste0(outdir, "/histogram.pdf"))
      hist(lens, 100, xlim = range(500,1600))
      dev.off()
      
      filts <- file.path(indir, "filtered", basename(fns))
      track <- filterAndTrim(fns, filts, minQ=minQ.a, minLen=1000, maxLen=1600, maxN=0, rm.phix=FALSE, maxEE=maxEE.a, trimLeft=trimLeft.a, truncQ=truncQ.a)
      write.table(track, paste0(outdir, "/filter_and_trim_stats.tsv"), sep="\t", quote=FALSE, col.names=NA)

      drp <- derepFastq(filts, verbose=TRUE)
      err <- learnErrors(drp, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=num_threads)
      saveRDS(err, file.path(outdir, "err.rds"))
      pdf(file=paste0(outdir, "/err.pdf"))
      plotErrors(err)
      dev.off()
      dd <- dada(drp, err=err, BAND_SIZE=32, multithread=TRUE)
      saveRDS(dd, file.path(outdir, "dd.rds"))
     
      st <- makeSequenceTable(dd)
      saveRDS(st, paste0(outdir, "/seqtab.rds"))
      st2 = data.frame(t(st), check.names=FALSE)
      curr_colnames = colnames(st2)
      curr_colnames = gsub("_R1.fastq.gz", "", curr_colnames)
      curr_colnames = gsub("_R1.fastq", "", curr_colnames)
      colnames(st2) = curr_colnames
      # Construct sequence table
      st2$otu_id = row.names(st2)
      st2 = st2[,c(ncol(st2), 1:(ncol(st2)-1))]
      colnames(st2)[1] = "#FEATURE_ID"
      write.table(st2, paste0(outdir, "/all.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
      
      # Scan for chimeras
      seqtab = NULL
      print("Scanning for chimeras and writing chimera-less table...")
      seqtab = makeSequenceTable(dd)
      seqtab = removeBimeraDenovo(seqtab, method="consensus", multithread=num_threads, minFoldParentOverAbundance=3.5)
      saveRDS(seqtab, paste0(outdir, "/seqtab_nochimera.rds"))
      
      seqtab =  data.frame(t(seqtab), check.names=FALSE)
      print(head(seqtab))
      seqtab$otu_id = row.names(seqtab)
      seqtab = seqtab[,c(ncol(seqtab), 1:(ncol(seqtab)-1))]
      colnames(seqtab)[1] = "#FEATURE_ID"
      write.table(seqtab, paste0(outdir, "/all_nochimera.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
   
   }else if(type == "spe"){
      outdir_filt = file.path(outdir, "filtered") 
      fastqFs <- sort(list.files(indir, pattern="_R1.fastq", full.names = FALSE))
      fastqRs <- sort(list.files(indir, pattern="_R2.fastq", full.names = FALSE))
      if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
      
      # Filter and trim
      print("executing filterAndTrim() function...")
      track = filterAndTrim(fwd=file.path(indir, fastqFs), filt=file.path(outdir_filt, fastqFs),
                    rev=file.path(indir, fastqRs), filt.rev=file.path(outdir_filt, fastqRs),
                    trimLeft=trimLeft.a, maxEE=maxEE.a, truncQ=truncQ.a, maxN=0, rm.phix=TRUE,
                    compress=FALSE, verbose=TRUE, multithread=TRUE, matchIDs=FALSE)
      
      write.table(track, paste0(outdir, "/filter_and_trim_stats.tsv"), sep="\t", quote=FALSE, col.names=NA)
      filtFs <- list.files(outdir_filt, pattern="_R1.fastq", full.names = TRUE)
      filtRs <- list.files(outdir_filt, pattern="_R2.fastq", full.names = TRUE)
      sample.namesF <- sapply(gsub("_R1.fastq", "", basename(filtFs)), `[`, 1)
      sample.namesR <- sapply(gsub("_R2.fastq", "", basename(filtRs)), `[`, 1)

      print("sample.namesF:")
      print(sample.namesF)
      print("sample.namesR:")
      print(sample.namesR)

      if(!identical(sample.namesF, sample.namesR)) stop("Forward and reverse files do not match.")
      names(filtFs) <- sample.namesF
      names(filtRs) <- sample.namesF
      sample.names = sample.namesF
      set.seed(100)
      # Learn forward error rates
      errF <- learnErrors(filtFs, nbases=1e8, multithread=num_threads, verbose=TRUE)
      # Learn reverse error rates
      errR <- learnErrors(filtRs, nbases=1e8, multithread=num_threads, verbose=TRUE)
      pdf(file=paste0(outdir, "/errF.pdf"))
      plotErrors(errF)
      dev.off()
      pdf(file=paste0(outdir, "/errR.pdf"))
      plotErrors(errR)
      dev.off()
      # Sample inference and merger of paired-end reads
      mergers <- vector("list", length(sample.names))
      names(mergers) <- sample.names
      for(sam in sample.names) {
         cat("Processing:", sam, "\n")
         cat("Derep...R1\n")
         derepF <- derepFastq(filtFs[[sam]])
         cat("dada...R1\n")
         ddF <- dada(derepF, err=errF, multithread=num_threads)
         cat("Derep...R2\n")
         derepR <- derepFastq(filtRs[[sam]])
         cat("dada...R2\n")
         ddR <- dada(derepR, err=errR, multithread=num_threads)
         cat("Merging pairs...\n")
         merger <- mergePairs(ddF, derepF, ddR, derepR, minOverlap=min_overlap, maxMismatch=1, verbose=TRUE)
         mergers[[sam]] <- merger
      }
      rm(derepF); rm(derepR)

      # Construct sequence table
      saveRDS(makeSequenceTable(mergers), paste0(outdir, "/seqtab.rds"))
      seqtab <- data.frame(t(makeSequenceTable(mergers)), check.names=FALSE)
      print(head(seqtab))
      seqtab$otu_id = row.names(seqtab)
      seqtab = seqtab[,c(ncol(seqtab), 1:(ncol(seqtab)-1))]
      colnames(seqtab)[1] = "#FEATURE_ID"
      write.table(seqtab, paste0(outdir, "/all.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

      # Scan for chimeras
      seqtab = NULL
      print("Scanning for chimeras and writing chimera-less table...")
      seqtab = makeSequenceTable(mergers)
      seqtab = removeBimeraDenovo(seqtab, method="consensus", multithread=num_threads)
      saveRDS(seqtab, paste0(outdir, "/seqtab_nochimera.rds"))
      
      seqtab =  data.frame(t(seqtab), check.names=FALSE)
      print(head(seqtab))
      seqtab$otu_id = row.names(seqtab)
      seqtab = seqtab[,c(ncol(seqtab), 1:(ncol(seqtab)-1))]
      colnames(seqtab)[1] = "#FEATURE_ID"
      write.table(seqtab, paste0(outdir, "/all_nochimera.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
      
   }else if(type == "sse"){
      outdir_filt = file.path(outdir, "filtered") 
      fastqFs <- sort(list.files(indir, pattern="_R1.fastq", full.names = FALSE))
      # Filter and trim
      track = filterAndTrim(fwd=file.path(indir, fastqFs), filt=file.path(outdir_filt, fastqFs),
                    trimLeft=trimLeft.a, maxEE=maxEE.a, truncQ=truncQ.a, maxN=0, rm.phix=TRUE,
                    compress=FALSE, verbose=TRUE, multithread=num_threads)
      write.table(track, paste0(outdir, "/filter_and_trim_stats.tsv"), sep="\t", quote=FALSE, col.names=NA)
      
      filtFs <- list.files(outdir_filt, pattern="_R1.fastq", full.names = TRUE)
      sample.namesF <- sapply(gsub("_R1.fastq", "", basename(filtFs)), `[`, 1)
      names(filtFs) <- sample.namesF
      set.seed(100)
      # Learn forward error rates
      errF <- learnErrors(filtFs, nbases=1e8, multithread=num_threads)
      pdf(file=paste0(outdir, "/errF.pdf"))
      plotErrors(errF)
      dev.off()
      # Sample inference and merger of paired-end reads
      mergers <- vector("list", length(sample.namesF))
      names(mergers) <- sample.namesF
      for(sam in sample.namesF) {
         cat("Processing:", sam, "\n")
         derepF <- derepFastq(filtFs[[sam]])
         ddF <- dada(derepF, err=errF, multithread=num_threads)
         mergers[[sam]] <- ddF
      }

      # Construct sequence table and remove chimeras
      saveRDS(makeSequenceTable(mergers), paste0(outdir, "/seqtab.rds"))
      seqtab <- data.frame(t(makeSequenceTable(mergers)), check.names=FALSE)
      seqtab$otu_id = row.names(seqtab)
      seqtab = seqtab[,c(ncol(seqtab), 1:(ncol(seqtab)-1))]
      colnames(seqtab)[1] = "#FEATURE_ID"
      write.table(seqtab, paste0(outdir, "/all.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
      saveRDS(seqtab, paste0(outdir, "/seqtab.rds"))
      
      # Scan for chimeras
      seqtab = NULL
      print("Scanning for chimeras and writing chimera-less table...")
      seqtab = makeSequenceTable(mergers)
      seqtab = removeBimeraDenovo(seqtab, method="consensus", multithread=num_threads)
      saveRDS(seqtab, paste0(outdir, "/seqtab_nochimera.rds"))
      
      seqtab =  data.frame(t(seqtab), check.names=FALSE)
      print(head(seqtab))
      seqtab$otu_id = row.names(seqtab)
      seqtab = seqtab[,c(ncol(seqtab), 1:(ncol(seqtab)-1))]
      colnames(seqtab)[1] = "#FEATURE_ID"
      write.table(seqtab, paste0(outdir, "/all_nochimera.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
   
   }else if(type == "sse_nc2"){
      outdir_filt = file.path(outdir, "filtered") 
      fastqRs <- sort(list.files(indir, pattern="_R2.fastq", full.names = FALSE))
      # Filter and trim
      track = filterAndTrim(rev=NULL, fwd=file.path(indir, fastqRs), filt=file.path(outdir_filt, fastqRs),
                    trimLeft=trimLeft.a, maxEE=maxEE.a, truncQ=truncQ.a, maxN=0, rm.phix=TRUE,
                    compress=FALSE, verbose=TRUE, multithread=num_threads)
      write.table(track, paste0(outdir, "/filter_and_trim_stats.tsv"), sep="\t", quote=FALSE, col.names=NA)
      
      filtFs <- list.files(outdir_filt, pattern="_R2.fastq", full.names = TRUE)
      sample.namesF <- sapply(gsub("_R2.fastq", "", basename(filtFs)), `[`, 1)
      names(filtFs) <- sample.namesF
      set.seed(100)
      # Learn forward error rates
      errF <- learnErrors(filtFs, nbases=1e8, multithread=num_threads)
      pdf(file=paste0(outdir, "/errF.pdf"))
      plotErrors(errF)
      dev.off()
      # Sample inference and merger of paired-end reads
      mergers <- vector("list", length(sample.namesF))
      names(mergers) <- sample.namesF
      for(sam in sample.namesF) {
         cat("Processing:", sam, "\n")
         derepF <- derepFastq(filtFs[[sam]])
         ddF <- dada(derepF, err=errF, multithread=num_threads)
         mergers[[sam]] <- ddF
      }

      # Construct sequence table and remove chimeras
      saveRDS(makeSequenceTable(mergers), paste0(outdir, "/seqtab.rds"))
      seqtab <- data.frame(t(makeSequenceTable(mergers)), check.names=FALSE)
      seqtab$otu_id = row.names(seqtab)
      seqtab = seqtab[,c(ncol(seqtab), 1:(ncol(seqtab)-1))]
      colnames(seqtab)[1] = "#FEATURE_ID"
      write.table(seqtab, paste0(outdir, "/all.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
      saveRDS(seqtab, paste0(outdir, "/seqtab.rds"))
      
      # Scan for chimeras
      seqtab = NULL
      print("Scanning for chimeras and writing chimera-less table...")
      seqtab = makeSequenceTable(mergers)
      seqtab = removeBimeraDenovo(seqtab, method="consensus", multithread=num_threads)
      saveRDS(seqtab, paste0(outdir, "/seqtab_nochimera.rds"))
      
      seqtab =  data.frame(t(seqtab), check.names=FALSE)
      print(head(seqtab))
      seqtab$otu_id = row.names(seqtab)
      seqtab = seqtab[,c(ncol(seqtab), 1:(ncol(seqtab)-1))]
      colnames(seqtab)[1] = "#FEATURE_ID"
      write.table(seqtab, paste0(outdir, "/all_nochimera.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
   }
} 

usage=function(errM) {
        cat("\nUsage : Rscript dada2.R [option] <Value>\n")
        cat("       -i        : indir\n")
        cat("       -o        : outdir\n")
        cat("       -t        : type\n")
        cat("       -p        : num_threads\n")
        cat("       -l        : trimLeft\n")
        cat("       -m        : maxEE\n")
        cat("       -c        : truncQ\n")
        cat("       -q        : minQ\n")
        cat("       -v        : min_overlap\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 9) {
	usage("missing arguments")
}
print("ARGs")
print(ARG)

## get arg variables
for (i in 1:length(ARG)) {
	if(ARG[i] == "-i") {
		indir =     ARG[i+1]
	}else if (ARG[i] == "-o") {
		outdir =    ARG[i+1]
	}else if( ARG[i] == "-t"){
	   type =       ARG[i+1]
	}else if( ARG[i] == "-p"){
      num_threads = ARG[i+1]
	}else if( ARG[i] == "-l"){
	   trimLeft.a = ARG[i+1]
	}else if( ARG[i] == "-m"){
	   maxEE.a =    ARG[i+1]
	}else if( ARG[i] == "-c"){
	   truncQ.a =   ARG[i+1]
	}else if( ARG[i] == "-q"){
	   minQ.a =     ARG[i+1]
	}else if( ARG[i] == "-v"){
	   min_overlap = ARG[i+1]
	}

}

doDada2(indir, outdir, type, num_threads, trimLeft.a, maxEE.a, truncQ.a, minQ.a, min_overlap)

