#!/usr/bin/env Rscript

# studies to focus on
studies <- commandArgs(trailingOnly=TRUE)
if(length(studies) > 1){
  print("OOPS! Enter one project only: e.g. SKCM. Or, do not pass any parameter to process non-SKCM and non-LAML data.")
}else{
  # load required package
  if (!require('data.table')) {install.packages("data.table")}
  if (!require('tidyverse')) {install.packages("tidyverse")}
  suppressMessages(library(data.table))
  suppressMessages(library(tidyverse))

  # samples by project categories
  df_gene_all <- fread("https://tau.cmmt.ubc.ca/cSurvival/project_data/977/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv.header",nrows = 0) %>% colnames(.)
  if(studies == "SKCM"){
    infiles <- "df_gene.all_good.skcm_0106dup"; pat0 <- "^01|06"
  }else if(studies == "LAML"){
    infiles <- "df_gene.all_good.laml_0309dup"; pat0 <- "^03|09"
  }else{
    infiles <- "df_gene.all_good.non_skcmlaml_01dup"; pat0 <- "^01"
  }
  
  # the indices of primary tumors
  df_gene_all_01 <- sapply(df_gene_all, function(x){if(grepl(pat0,strsplit(x,"-")[[1]][4])){return(x)}else{return(NULL)}})
  df_gene_all_01 <- Filter(Negate(is.null), df_gene_all_01)
  
  if(studies != "LAML"){
    # extract duplicated samples
    df_gene_dup <- fread(paste0("https://tau.cmmt.ubc.ca/cSurvival/project_data/977/",infiles),header = F) %>% .[["V2"]]
    
    # the pattern to search the indices of duplicated samples
    pat <- paste0(df_gene_dup,collapse="|")
    pat <- (grepl(paste0("^",pat),df_gene_all) & (df_gene_all %in% df_gene_all_01))
  }
  
  if(studies != "LAML"){
    # the output files
    ofile <- paste0("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.good.",gsub("^df_gene.all_good.","",infiles),".tsv")
    
    # extract duplicated samples' data and calculate geometric means
    data <- fread("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",select = which(pat))
    data <- log10(data + 1)
    odd_indexes<-seq(1,ncol(data),2)
    for(i in odd_indexes){
      data[[i]] <- apply(data[,i:(i+1)],1,mean)
    }
    data <- dplyr::select(data,all_of(odd_indexes))
    data <- 10^data - 1
    data <- t(data)
    patient_ids <-sapply(rownames(data),function(x) paste0(strsplit(x,"-")[[1]][1:3],collapse = "-"))
    rownames(data) <- patient_ids
    fwrite(as.data.frame(data),ofile,sep="\t",row.names = T, col.names = F,na = NA,quote=F)
  }
}