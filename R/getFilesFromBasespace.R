#!/usr/bin/env Rscript

# Function that fetches files from basespace given an input project id.
# Author: Julien Tremblay - jtremblay514@gmail.com	
getFilesFromBasespace <- function(PROJECT_ID){	
	library(BaseSpaceR)
	
   ### OLD access token for julie.champagne@cnrc.nrc.gc.ca
   #ACCESS_TOKEN <- '09876ee5b7d24dfbb95dbcb138667277'
   
   ### Current access token for sylvie.sanschagrin@nrc-cnrc.gc.ca
   ACCESS_TOKEN <- '2e560b542f284938add49077bc0fb9b6'
	#ID=8782792
	#TOKEN=09876ee5b7d24dfbb95dbcb138667277
	#PROJECT_ID <- '18733732'  ## Get proj ID from url of the project

   PROJECT_ID = as.character(PROJECT_ID)	
	aAuth <- AppAuth(access_token = ACCESS_TOKEN)
	selProj <- Projects(aAuth, id = PROJECT_ID, simplify = TRUE) 
	sampl <- listSamples(selProj, limit= 1000)
	inSample <- Samples(aAuth, id = Id(sampl), simplify = TRUE)
	for(s in inSample){ 
	   f <- listFiles(s, Extensions = ".gz")
	   #print(Name(f))
	   #file_name = paste0("./Data/Intensities/BaseCalls/", Name(f)) 
	
	   for(i in 1:length(f)){
	      curr_file_name = paste0("./Data/Intensities/BaseCalls/", Name(f[i]))
	      if(file.exists(curr_file_name) && (file.info(curr_file_name)$size[1] > 0)){ 
	            print(paste0("File: ", curr_file_name, " already exists and is not empty..."))
	      }else{
	         getFiles(aAuth, id= Id(f[i]), destDir = './', verbose = TRUE)
	         print(paste0("Downloading File: ", curr_file_name, " ..."))
	      }
	   }
	}
}

usage=function(errM) {
        cat("\nUsage : Rscript getFilesFromBasespace.R [option] <Value>\n")
        cat("       -p        : project_id\n")
}

ARG = commandArgs(trailingOnly = T)

if(length(ARG) < 1) {
    usage("missing arguments")
}

## get arg variables
for (i in 1:length(ARG)) {
    if (ARG[i] == "-p") {
        project_id=ARG[i+1]
    }
}

getFilesFromBasespace(project_id)
