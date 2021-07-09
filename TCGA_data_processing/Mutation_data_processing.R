#!/usr/bin/env Rscript

# studies to focus on
studies <- commandArgs(trailingOnly=TRUE)

# load required package
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require('org.Hs.eg.db')) {BiocManager::install("org.Hs.eg.db")}
if (!require('data.table')) {install.packages("data.table")}
if (!require('tidyverse')) {install.packages("tidyverse")}
library(org.Hs.eg.db)
library(data.table)
library(tidyverse)

# avaialbel TCGA studies
tcga_projects <- fread("https://tau.cmmt.ubc.ca/cSurvival/project_data/977/projects",header=F)
tcga_projects <- tcga_projects[["V1"]]

# read in MAF mutations, extract essential data
infiles <- "mc3.v0.2.8.PUBLIC.good_01.maf"
if("SKCM" %in% studies){
  infiles <- c(infiles,"mc3.v0.2.8.PUBLIC.good_06.maf")
}
if("LAML" %in% studies){
  infiles <- c(infiles,"mc3.v0.2.8.PUBLIC.good_0309.maf")
}
l <- lapply(infiles,function(y){
  fread(y,select = c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode"))
})
muts <- rbindlist(l, use.names = T)

# convert sample bar codes to patient IDs
muts[["Tumor_Sample_Barcode"]] <-
  sapply(muts[["Tumor_Sample_Barcode"]], function(x){
    x <- strsplit(x, "-")[[1]][1:3]
    paste0(x, collapse = "-")
  })

# convert gene IDs
egSYMBOL <- toTable(org.Hs.egSYMBOL)
muts <- muts %>% dplyr::left_join(egSYMBOL, by=c("Hugo_Symbol"="symbol"))
muts[["Hugo_Symbol"]] <- paste0(muts[["Hugo_Symbol"]],"|",muts[["gene_id"]])
muts <- muts %>% dplyr::select(-gene_id)

# extract data per study
for(i in tcga_projects){
  if(length(studies)==0){
    if(i == "SKCM" | i == "LAML"){next}
  }else{
    if(!(i %in% studies)){next}
  }
  odir <- paste0("./TCGA-",i)
  if(!dir.exists(odir)){dir.create(odir)}
  ofile <- paste0(odir,"/df_snv_class_977.csv")
  i_ids <- fread(paste0(i,".testt"),header=F) %>% .[["V1"]]
  i_muts <- muts %>% dplyr::filter(Tumor_Sample_Barcode %in% i_ids)
  muts_genes <- c("patient_id",unique(i_muts[["Hugo_Symbol"]]))
  fwrite(as.list(muts_genes %>% paste0(collapse = ",")), ofile, quote=F)
  muts_patients <- unique(i_muts[["Tumor_Sample_Barcode"]])
  for(patient in muts_patients){
    e_muts <- i_muts %>% dplyr::filter(Tumor_Sample_Barcode == patient)
    genes <- unique(e_muts[["Hugo_Symbol"]])
    e_mutss <- lapply(genes, function(x){
      e_muts %>% dplyr::filter(Hugo_Symbol == x) %>%
        .[["Variant_Classification"]] %>% unique() %>% paste0(collapse = "|")
    })
    names(e_mutss) <- genes
    w_muts <- lapply(muts_genes[-1], function(x){
      if(x %in% names(e_mutss)){
        e_mutss[[x]]
      }else{
        ""
      }
    }) %>% paste0(collapse = ",")
    w_muts <- paste0(patient,",",w_muts)
    fwrite(as.list(w_muts), ofile, quote=F, append = T)
  }
}