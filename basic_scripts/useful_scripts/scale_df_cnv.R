#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Please enter the input file prefix (e.g. df_mir).", call.=FALSE)
}
library(data.table)

infile <- args[1]

test_path <- paste0(infile,".csv")

scaled_df <- function(df_path,sep=","){
  genes <- fread(df_path, sep=sep, nrows = 0)
  genes <- colnames(genes)[-1]
  df <- fread(df_path, sep=sep, select = "patient_id")
  
  ofile <- paste0(infile,"_scale.csv")
  unlink(ofile)
  ogenes <- "patient_id"
  
  for(i in seq_along(genes)){
    gene <- genes[i]
    data <- fread(df_path, sep=sep, select = c("patient_id",gene))
    data <- data[[gene]]
    yn <- table(data[!is.na(data)])
    if(length(yn) == 1){
      
    }else if((length(yn) == 2 & min(yn)>10) | (length(yn) > 2 & sort(yn)[[2]]>10)){
      ogenes <- c(ogenes,gene)
    }
  }
  fwrite(as.list(paste0(ogenes,collapse=",")),ofile, sep = ",", quote = F)
}

scaled_df(test_path)
