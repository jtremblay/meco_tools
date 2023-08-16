#!/usr/bin/env Rscript

# Function solving polytope equation for Ray assembler input.
solvePolytopeDegreeBruteForce <- function(ncores,top){
  r = expand.grid(c(1:10000),c(1:1000))
  r$v = r[,1]^r[,2]
  colnames(r) = c('alphabetSize','wordLength','vertices')
  r$degree = (r[,1]-1) * r[,2]
  r$diff =  as.numeric(ncores) - r$vertices
  r = r[r$diff>=0 & r$wordLength>1,]
  r = r[order(r$diff),]
  return(head(r,top))
}

usage=function(errM) {
        cat("\nUsage : Rscript findPolytope.R [option] <Value>\n")
        cat("       -c        : Number of cores\n")
        cat("       -n        : Number of nodes\n")
        cat("       -o        : outfile (csv)\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 4) {
	usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
	if(ARG[i] == "-o") {
		outfile=ARG[i+1]
	}else if (ARG[i] == "-c") {
		cores=ARG[i+1]
	}else if (ARG[i] == "-n") {
		nodes=ARG[i+1]
	}
}

#makeplots(data, outdir, prefix)
values = solvePolytopeDegreeBruteForce(cores, nodes)
write.csv(values, outfile)

