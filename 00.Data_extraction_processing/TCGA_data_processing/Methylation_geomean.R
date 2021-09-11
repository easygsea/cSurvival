#!/usr/bin/env Rscript

dir_path_o <- getwd()
dir_path <- file.path(dir_path_o, "methylation") 
setwd(dir_path)
met_files <- list.files(pattern = "^TCGA-.*tsv$")

geo_mean <- function(x){if(!all(is.na(x))){x <- log(x)};mean(x,na.rm=T)}
for(i in seq_along(met_files)){
  file <- met_files[i]
  data <- fread(file, sep="\t")
  sample <- unlist(data[1,1])
  sample <- sapply(sample, function(x) paste0(strsplit(x,"-")[[1]][1:3],collapse = "-")) %>% unname()
  data <- apply(data[,-1], 2, geo_mean)
  ofile <- paste0(file,".geomean")
  fwrite(as.list(paste0(c(sample,data),collapse = "\t")),ofile)
}

setwd(dir_path_o)