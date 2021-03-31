library(tidyverse)
library(RCurl) # to download the specific folder
library(tools) # to check the file name extensions
library(R.utils) # unzip the .gz files
library(survival) # to do the survival analysis
library(survminer) # to plot the survival analysis nicer
library(data.table)
library(TCGAbiolinks)
library(microbenchmark) # test the time spent on a command
library(fgsea) # use the function gmtpathways
library("org.Hs.eg.db") # generate the id conversion table
library(maftools)

select <- dplyr::select
# setwd("C:/Users/15067/Desktop/WORK!!!!!!!!1/gdc_data")
# setwd("~/Desktop/cSurvival/basic_scripts")
# setwd("/Users/jeancheng/Documents/KMplot/LUAD_V2_TCGAbiolink")
df_gdc_projects <- getGDCprojects() %>%
  filter(grepl("TCGA", id))


# --------------------

# generate the table containing the full name of all genes using the id conversion table
# input: a vector of gene id, each element starts with "ENSG..."
# output: a full name conversion table
generate_full_name_table <- function(gene_ids_ensg){
  # create individual tables using org.Hs
  egENS <- toTable(org.Hs.egENSEMBL)
  egSYMBOL <- toTable(org.Hs.egSYMBOL)
  
  # bind the tables
  id_table <- egENS %>% left_join(egSYMBOL, by = "gene_id")
  
  # extract the gene ids from df_gene and find their names from id conversion table
  gene_ids <- as_tibble_col(gene_ids_ensg, column_name = "ensembl_id")
  gene_ids_table <- left_join(gene_ids, id_table, by="ensembl_id")
  
  # extract the ids that do not exist in the id conversion table
  gene_na_table <- gene_ids_table %>%
    filter(is.na(gene_id)) %>%
    mutate(full_name = ensembl_id)
  # extract the ids that exist in the id conversion table
  gene_not_na <- gene_ids_table %>%
    filter(!is.na(gene_id))
  # summarise the duplicated gene_ids and symbols in the df
  unique_ids_symbols <-
    aggregate(x= tibble(gene_not_na["gene_id"], gene_not_na["symbol"]),
              by=gene_not_na["ensembl_id"],
              FUN = function(X) paste(unique(X), collapse="|")) %>%
    mutate(full_name = paste(symbol, gene_id, ensembl_id, sep = "|"))
  # the table that contains all the full name of "ENSG...."
  full_name_table <- rbind(unique_ids_symbols, gene_na_table)
}
# the name of the cancer project
project_name <- "TCGA-LUAD"
# the query to download the data related to .cnv file
query_cnv <- try(
  GDCquery(project = project_name,
           data.category = "Copy Number Variation",
           data.type = "Gene Level Copy Number Scores",
           experimental.strategy = "Genotyping Array",
           platform = "Affymetrix SNP 6.0")
)
# download the corresponding cnv data
GDCdownload(query_cnv, method = "api", directory = paste0(project_name, "/", "cnv_data"))
# read the file and clean up the gene ids
filename_cnv <- list.files(pattern = query_cnv[[1]][[1]]$file_name[1], recursive = T)
df_cnv <- fread(filename_cnv) %>%
  mutate(`Gene Symbol` = gsub("[.]\\d+$","",`Gene Symbol`)) %>%
  rename(`Gene Symbol` = 'full_name')
# convert the gene ids into the full names of genes
full_names_cnv <- left_join(as_tibble_col(df_cnv$full_name, column_name = "full_name"),
                                 generate_full_name_table(df_cnv$full_name),
                                 by = c("full_name" = "ensembl_id"))$full_name.y
df_cnv <- df_cnv %>%
  mutate(full_name = full_names_cnv)

#------------------------------------------------------------------------------------------------





# create a query that downloads all the snv data
query_snv <- try(
                  GDCquery(project = project_name,
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      experimental.strategy = "WXS")
                )
if(inherits(query_snv, "try-error")){next}
# download them
GDCdownload(query_snv, method = "api", directory = paste0(project_name, "/", "snv_data"))
# find the exact file that we want
special_names <- c("mutect", "varscan", "somaticsniper", "muse")
for(name_special in special_names){
    # name_special <- "mutect"
    # the path to output the data frames related to snv
    df_snv_type_path = paste0(project_name, "/df_snv_type_", name_special, ".csv")
    df_snv_class_path = paste0(project_name, "/df_snv_class_", name_special, ".csv")
    # read the maf
    filename_snv <- list.files(pattern = grep(name_special, query_snv[[1]][[1]]$file_name,value = T), recursive = T)[1]
    df_maf <- try(read.maf(maf  = filename_snv))
    # extract the useful information and generate the column names
    df_snv <- try(
      df_maf@data[,c("Hugo_Symbol", "Entrez_Gene_Id", "Gene", 
                             "Tumor_Sample_Barcode","Variant_Type", "Variant_Classification")] %>%
      mutate(full_name = paste(Hugo_Symbol, Entrez_Gene_Id, Gene, sep = "|")) %>%
      mutate(Tumor_Sample_Barcode = str_sub(Tumor_Sample_Barcode, start = 1L, end = 12L))
    )
    # if the data is read successfully, delete the older version of csv
    if(!inherits(df_snv, "try-error")){
      unlink(df_snv_type_path, recursive = T)
      unlink(df_snv_class_path, recursive = T)
    } else {next}
      
    # write the first line of the csvs
    write(paste0("patient,", paste(unique(df_snv$full_name), collapse = ",")), file = df_snv_type_path, append = T)
    write(paste0("patient,", paste(unique(df_snv$full_name), collapse = ",")), file = df_snv_class_path, append = T)
    
    # the list of data frames that contain the infomation of each patient
    snv_dfs <- split(df_snv, df_snv$Tumor_Sample_Barcode)
    
    for(i in seq_along(snv_dfs)){
      patient_snv <-left_join(tibble(full_name = unique(df_snv$full_name)), snv_dfs[[i]][!duplicated(snv_dfs[[i]]$full_name), ], by = "full_name")
      # write the information to corresponding csvs
      write(paste0(names(snv_dfs)[i], ",", paste(patient_snv$Variant_Type, collapse = ",")), file = df_snv_type_path, append = T)
      write(paste0(names(snv_dfs)[i], ",", paste(patient_snv$Variant_Classification, collapse = ",")), file = df_snv_class_path, append = T)
    }

}


df_blablabla <- fread(df_snv_type_path)




project_name <- "TCGA-COAD"
# the path to write all the csvs
df_survival_path <- paste0(project_name, "/df_survival.csv")
df_gene_path = paste0(project_name, "/", "df_gene.csv")
df_gene_scale_path = paste0(project_name, "/", "df_gene_scale.csv")
# The path for reading and writing the patient data frames
dir.create(project_name)
# the query to download the gene information of all patients
query_gene_counts <- try(
                      GDCquery(project = project_name,
                              data.category = "Transcriptome Profiling",
                              experimental.strategy = "RNA-Seq",
                              data.type = "Gene Expression Quantification",
                              workflow.type = "HTSeq - FPKM",
                              access = "open")
                      )
if(!inherits(query_gene_counts, "try-error")){
  unlink(paste0(project_name, "/", "gene_data"), recursive = T)
} else {next}
GDCdownload(query_gene_counts, method = "api", directory = paste0(project_name, "/", "gene_data"))

# the query to download the survival information of all patients
query_patient_clinical <- try(
                            GDCquery(project = project_name,
                                   data.category = "Clinical",
                                   data.type = "Clinical Supplement",
                                   data.format = "BCR Biotab")
                          )
if(!inherits(query_patient_clinical, "try-error")){
  unlink(paste0(project_name, "/", "patient_data"), recursive = T)
} else {next}
GDCdownload(query_patient_clinical, method = "api", directory = paste0(project_name, "/", "patient_data"))


filename_patient <- list.files(pattern = grep("clinical_patient", query_patient_clinical[[1]][[1]]$file_name,value = T), recursive = T)[1]

cases_name_complete <- query_gene_counts[[1]][[1]]$cases
# the vector contains the shorten name of all patient ids
cases_name_shorten <-
  map_chr(cases_name_complete, function(x){
    str_sub(x, start = 1L, end = 12L)
  })
# the list to store the duplicated and non-duplicated patient information
cases_info <- list()
for(i in seq_along(cases_name_shorten)){
  patient_name <- cases_name_shorten[i]
  file_name <- query_gene_counts[[1]][[1]]$file_name[i]
  if(is.null(cases_info[[patient_name]])){
    # create a new list of unexisted patient
    cases_info[[patient_name]] <- list(file_name)
  } else {
    # store the filename in the list
    cases_info[[patient_name]][[length(cases_info[[patient_name]])+1]] <- file_name
  }
}


# write the gene names into a file a file as the first line
filename_first_FPKM <- list.files(pattern = cases_info[[1]][[1]], recursive = TRUE)[1]
# sample_FPKM_name <- read_tsv(gzfile(filename_first_FPKM), col_names = FALSE)$X1
sample_FPKM_name <- try(read_tsv(gzfile(filename_first_FPKM), col_names = FALSE)$X1)
if(!inherits(sample_FPKM_name, "try-error")){
  unlink(df_gene_path, recursive = T)
} else {next}

write(paste0("patient_id,", paste(gsub("[.]\\d+$","",sample_FPKM_name), collapse = ",")), file = df_gene_path, append = T)
# go through all files and manipulate the duplicates
for(name in names(cases_info)){
  if(length(cases_info[[name]]) == 1){
    # when no duplicate samples of the patient, write a line to the output file
    filename_gene_info = list.files(pattern = cases_info[[name]][[1]], recursive = T)
    df_gene_write <- read_tsv(gzfile(filename_gene_info[1]), col_names = FALSE)
    df_gene_write <- left_join(tibble(X1 = sample_FPKM_name), df_gene_write, by = "X1") %>%
      mutate(X2 = log(as.numeric(X2)+1))
    write(paste0(name, ",", paste(df_gene_write$X2, collapse = ",")), file = df_gene_path, append = T)
    
  } else {
    # the list to store data frames with duplicated patient ids
    df_duplicated_list <- list()
    for(i in seq_along(cases_info[[name]])){
      df_duplicated_list[[i]] <- 
        read_tsv(gzfile(list.files(pattern = cases_info[[name]][[i]], recursive = T)[1]), col_names = FALSE) %>%
        mutate(X2 = log(as.numeric(X2)+1))
      names(df_duplicated_list[[i]])[names(df_duplicated_list[[i]]) == "X2"] <- paste0(name, "_", i)
      # print(list.files(pattern = cases_info[[name]][[i]], recursive = T)[1])
    }
    # calculate the means of all duplicates
    df_duplicated <- left_join(tibble(X1 = sample_FPKM_name),
                               Reduce(function(...) base::merge(..., by = "X1", all.x = T), df_duplicated_list),
                               by = "X1") %>%
                      mutate(mean = rowMeans(select(.,-X1),na.rm = T))
    write(paste0(name, ",", paste(df_duplicated$mean, collapse = ",")), file = df_gene_path, append = T)
    
    
  }
}

# create the df_survival
df_survival <- try(
  read_tsv(filename_patient, skip = 1, na = c("", "NA", "[Not Available]", "[Not Applicable]", "[Discrepancy]")) %>%
  select(patient_id = `bcr_patient_barcode`, death_days = `days_to_death`, followup_days = `days_to_last_followup`, diagnosis_days = `days_to_initial_pathologic_diagnosis`) %>%
  mutate(censoring_status = ifelse(is.na(death_days), 0, 1)) %>% # The status indicator, normally 0=alive, 1=dead. 
  filter(diagnosis_days == 0) %>%
  mutate(survival_days = as.numeric(ifelse(is.na(death_days), followup_days, death_days))) %>%
  filter(!is.na(survival_days))
)
if(!inherits(df_survival, "try-error")){
  unlink(df_survival_path, recursive = T)
} else {next}

# output the df_survival
fwrite(df_survival, file = df_survival_path)


# df_gene_scale <- read_csv(df_gene_path)
df_gene <- data.table::fread(df_gene_path)

# df_gene_scale <- data.table::fread(df_gene_path)
# df_gene <- df_gene_scale %>%
#   data.table::transpose(keep.names = "gene", make.names = "patient_id") %>% 
#   mutate_at(vars(!gene), as.numeric)
patient_ids <- df_gene$patient_id

# create individual tables using org.Hs
egENS <- toTable(org.Hs.egENSEMBL)
egSYMBOL <- toTable(org.Hs.egSYMBOL)

# bind the tables
id_table <- egENS %>% left_join(egSYMBOL, by = "gene_id")

# extract the gene ids from df_gene and find their names from id conversion table
gene_ids <- as_tibble_col(colnames(df_gene)[-1], column_name = "ensembl_id")
gene_ids_table <- left_join(gene_ids, id_table, by="ensembl_id")

# extract the ids that do not exist in the id conversion table
gene_na_table <- gene_ids_table %>%
  filter(is.na(gene_id)) %>%
  mutate(full_name = ensembl_id)
# extract the ids that exist in the id conversion table
gene_not_na <- gene_ids_table %>%
  filter(!is.na(gene_id))
# summarise the duplicated gene_ids and symbols in the df
unique_ids_symbols <- 
  aggregate(x= tibble(gene_not_na["gene_id"], gene_not_na["symbol"]), 
            by=gene_not_na["ensembl_id"], 
            FUN = function(X) paste(unique(X), collapse="|")) %>%
  mutate(full_name = paste(symbol, gene_id, ensembl_id, sep = "|"))

full_name_table <- rbind(unique_ids_symbols, gene_na_table)

# unique_symbols <- 
#   aggregate(x= gene_ids_table["gene_id"], 
#             by=gene_ids_table["ensembl_id"], 
#             FUN = function(X) length(unique(X)))
# filter(unique_symbols, gene_id == 7)
# max(unique_symbols$gene_id)

# scale the df_gene and generate df_gene_scale
df_gene_scale <- apply(df_gene[,-1], 2, scale)
# df_gene_scale<- df_gene_scale[,complete.cases(t(df_gene_scale))] 
# function to remove all na per column
not_all_na <- function(x) any(!is.na(x))
df_gene_scale <- try(
  as_tibble(df_gene_scale) %>%
  dplyr::select(where(not_all_na)) %>%
  mutate(patient_id = patient_ids) %>%
  dplyr::select(patient_id, everything())
)

colnames(df_gene_scale)[-1] <- 
  left_join(as_tibble_col(colnames(df_gene_scale)[-1], column_name = "ensembl_id"),
            full_name_table,
            by = "ensembl_id")$full_name

if(!inherits(df_gene_scale, "try-error")){
  unlink(df_gene_scale_path, recursive = T)
} else {next}

fwrite(df_gene_scale, file = df_gene_scale_path)
# fwrite(df_gene_scale_2, file = "df_gene_scale_2.csv")

# read into the survival data frame
df_survival <- fread(df_survival_path) %>% 
  mutate_at(vars(!patient_id), as.numeric)


# the function to get the data frame with lowest p.value ------------------
# gene_name: the name of a gene which exists in the gene data frame
# df_survival: the data frame having columns patient_id, survival_days, and censoring status
# df_gene: the gene df contains columns patient_id and gene(count of gene)
# output: a data frame containing columns survival_days, censoring_status and level, or
# the p.value, or the quantile, or the cutoff of gene level of the data frame
# attribute info = c("df", "pval", "quantile", "cutoff")
get_info_most_significant <- function(gene_name, rename=NULL, df_survival_previous = df_survival, df_gene_overall = df_gene_scale, min = .1, max = .9, step = .01){  
  # initiate quantiles according to margin and step values
  quantile_s = seq(min, max, by = step)
  # # generate the survival data frames with levels of gene percentages
  # df_survival_with_gene <- df_survival_previous %>%
  #   left_join(select(df_gene_overall, patient_id, gene = gene_name), by = "patient_id") %>%
  #   filter(!is.na(gene)) %>%
  #   mutate(gene = as.numeric(gene))
  # initialize the most significant p value and model
  least_p_value <- 1
  # df_most_significant <- tibble(df_survival_with_gene, level = 1)
  quantile_most_significant <- NULL

  # extract gene's data
  df_gene_s <- df_gene_overall %>% dplyr::select(patient_id, gene_name)
  patient_ids <- df_gene_s[["patient_id"]]
  
  # the quantiles we will use to define the level of gene percentages
  quantiles <- quantile(df_gene_s[[gene_name]], quantile_s)

  for(i in seq_along(quantiles)){
    # generate the data from for model first
    gene_quantiles <- df_gene_s[[gene_name]] %>% 
      lapply(function(x) ifelse(x > quantiles[[i]], "high", "low"))
    names(gene_quantiles) <- patient_ids
    # print(quantiles[[i]])
    # generate survival analysis df
    df_survival_quantile <- df_survival_previous
    df_survival_quantile$level <- gene_quantiles[match(df_survival_quantile$patient_id,names(gene_quantiles))]
    df_survival_quantile <- df_survival_quantile %>%
      dplyr::filter(level != "NULL")
    df_survival_quantile$level <- base::unlist(df_survival_quantile$level, recursive = F)
    # test if there is significant difference between high and low level genes
    surv_diff <- survdiff(Surv(survival_days, censoring_status) ~ level, data = df_survival_quantile)
    p_diff <- 1- pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
    if(p_diff < least_p_value){
      least_p_value = p_diff
      df_most_significant <- df_survival_quantile
      quantile_most_significant <- names(quantiles)[i]
      cutoff_most_significant <- quantiles[i]
    }
  }
  # if(info == "df"){
  #   return(df_most_significant)
  # } else if(info == "pval"){
  #   return(least_p_value)
  # } else if(info == "quantile"){
  #   return(quantile_most_significant)
  # } else if(info == "cutoff"){
  #   return(cutoff_most_significant)
  # } else{
  #   return("Your 'info' attribute must be values among df,pval,quantile,cutoff.")
  # }
  results <- list(list(
    df = df_most_significant,
    pval = least_p_value,
    quantile = quantile_most_significant,
    cutoff = cutoff_most_significant
  ))
  if(!is.null(rename)){
    names(results) <- rename
    
  }else{
    names(results) <- gene_name
    
  }
  return(results)
}

# get the level name consists of two genes in df_two_genes
get_level_name <- function(gene_name_1, gene_name_2){
  paste0(gene_name_1, "_", gene_name_2)
}

# create a data frame having the level combining two genes
# input: two gene names
get_df_two_genes <- function(gene_name_1, gene_name_2){
  df_two_genes <- get_info_most_significant(gene_name_1)[[1]][["df"]] %>%
    inner_join(dplyr::select(get_info_most_significant(gene_name_2)[[1]][["df"]], patient_id, level), by = "patient_id") %>%
    mutate("x_y_level" = paste0(`level.x`, "_", `level.y`))
  # df_two_genes[[get_level_name(gene_name_1, gene_name_2)]] <- paste0(df_two_genes$`level.x`, "_", df_two_genes$`level.y`)
  df_two_genes
}
gene_1 <- "ENSG00000198242"
results <- get_info_most_significant(gene_1)

# df_most_significant <- get_info_most_significant(gene_1, info = "df")
df_most_significant <- results[[1]][["df"]]
head(df_most_significant)
# (quantile_most_significant <- get_info_most_significant(gene_1, info = "quantile"))
# (least_p_value <- get_info_most_significant(gene_1, info = "pval"))
# (cutoff_most_significant <- get_info_most_significant(gene_1, info = "cutoff"))
km_fit <- survfit(Surv(survival_days, censoring_status) ~ level, data = df_most_significant) # type = "right"  
ggsurvplot(
  fit = km_fit,
  xlab = "Days",
  ylab = "Survival probability"
)

gene_2 <- "ENSG00000263089"
df_two_genes <- get_df_two_genes(gene_1, gene_2)
# the name of the level of the diagram
get_level_name(gene_1, gene_2)

km_fit_final <- survfit(Surv(survival_days, censoring_status) ~ x_y_level, data = df_two_genes)
head(df_two_genes)
sur_plot <- 
  ggsurvplot(
    fit = km_fit_final,
    xlab = "Days",
    ylab = "Survival probability",
    legend = "left",
    legend.title = get_level_name(gene_1, gene_2) # the legend title can be changed
    # palette = c("blue", "red", "green", "yellow")
    
  )
sur_plot
####---------------------------- TO DO later-------------------------
look_for_gene <- function(gene_name, df_all_genes){
  grep(gene_name, df_all_genes, value = TRUE)
  # grep(paste0("^",gene_name,"[|]"), df_all_genes, value = TRUE)
}
# the function to generate the data frame for survival analysis
# geneset: the name of a geneset; df_gmt_content: the df contains all gmtpathways; df_gene_scaled: the z-score df
generate_gene_gmt_df <- function(geneset, df_gmt_content = df_gmt, df_gene_scaled = df_gene_scale){
  gene_full_name <- c()
  for(i in seq_along(df_gmt_content[[geneset]])){
    gene <- look_for_gene(df_gmt_content[[geneset]][i], df_all_genes = colnames(df_gene_scaled)[-1])
    if(length(gene) > 0){
      gene_full_name[i] <- gene
    }
  }
  gene_full_name <- gene_full_name[!is.na(gene_full_name)]
  
  tibble(gene_full_name)
  df_gene_gmt <- df_gene_scaled %>%
    select(c("patient_id", all_of(gene_full_name))) %>%
    mutate(signature_score = rowMeans(select(.,all_of(gene_full_name)), na.rm = T)) %>%
    select(patient_id, signature_score, everything())
  
  patient_ids <- df_gene_gmt$patient_id
  gene_list <- colnames(df_gene_gmt)[-1]
  df_gene_gmt <- t(df_gene_gmt)[-1,] %>% as_tibble() %>% mutate_at(vars(everything()), as.numeric)
  colnames(df_gene_gmt) <- patient_ids
  df_gene_gmt$gene <- gene_list
  df_gene_gmt <- df_gene_gmt %>% dplyr::select(gene, everything())
  
  return(df_gene_gmt)
}

# # df_gene_example <- data.table::transpose(df_gene, keep.names = "gene", make.names = "patient_id")
# 
# # generate the data frame having z-scores as its elements
# gene_list <- df_gene$gene
# patient_ids <- colnames(df_gene)[-1]
# df_gene_scale <- apply(df_gene[,-1],1,scale) %>% t()
# rownames(df_gene_scale)<- gene_list
# df_gene_scale<- df_gene_scale[complete.cases(df_gene_scale),] 
# gene_list_new <- rownames(df_gene_scale)
# 
# df_gene_scale <- t(df_gene_scale) %>% as_tibble()
# df_gene_scale$patient_id <- patient_ids
# 
# df_gene_scale <- df_gene_scale %>% dplyr::select(patient_id, everything())

path_gmts <- list.files("./gmts")
df_gmt <- lapply(paste0("gmts/",path_gmts), function(x) gmtPathways(x)) %>% unlist(., recursive = F)
df_gmt[[1]]
####--------------------------This section is for TARGET ----------------------#
df_gdc_projects <- getGDCprojects() %>%
  filter(grepl("TARGET", id))

project_name <- "TARGET-ALL-P2"

# The path for reading and writing the patient data frames
dir.create(project_name)
# the path to output the data frames related to snv
df_snv_type_path = paste0(project_name, "/df_snv_type.csv")
df_snv_class_path = paste0(project_name, "/df_snv_class.csv")

# create a query that downloads all the snv data
query_snv <- try(
  GDCquery(project = project_name,
           data.category = "Simple Nucleotide Variation",
           data.type = "Masked Somatic Mutation",
           experimental.strategy = "WXS")
)
if(!inherits(query_snv, "try-error")){
  unlink(paste0(project_name, "/", "snv_data"), recursive = T)
  }else{next}
# download them
GDCdownload(query_snv, method = "api", directory = paste0(project_name, "/", "snv_data"))
# create a list that stores all the .mafs
df_target_snv_list <- list()
for(j in seq_along(query_snv[[1]][[1]]$file_name)){
  filename_target_snv <- list.files(pattern = query_snv[[1]][[1]]$file_name[j], recursive = T)
  df_maf_target_snv <- try(
    read.maf(maf = filename_target_snv
             ,useAll = F
             ,vc_nonSyn = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins"
                            ,"Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"
                            ,"Splice_Site","Translation_Start_Site"
                            ,"Silent","Splice_Region","Intron","5'UTR","RNA","3'UTR"
                            ,"5'Flank","3'Flank","IGR")
             )@data[,c("Hugo_Symbol", "Entrez_Gene_Id", "Gene", 
                                          "Tumor_Sample_Barcode","Variant_Type", "Variant_Classification")] %>%
      mutate(full_name = paste(Hugo_Symbol, Entrez_Gene_Id, Gene, sep = "|")) %>%
      mutate(Tumor_Sample_Barcode = str_sub(Tumor_Sample_Barcode, start = 1L, end = 16L))
  )
  if(!inherits(df_maf_target_snv, "try-error")){
    df_target_snv_list[[j]] <- df_maf_target_snv
  }else{next}
}
# combine together all the .mafs
df_snv <- try(bind_rows(df_target_snv_list))

# if the data is read successfully, delete the older version of csv
if(!inherits(df_snv, "try-error")){
  unlink(df_snv_type_path, recursive = T)
  unlink(df_snv_class_path, recursive = T)
} else {next}

# write the first line of the csvs
write(paste0("patient,", paste(unique(df_snv$full_name), collapse = ",")), file = df_snv_type_path, append = T)
write(paste0("patient,", paste(unique(df_snv$full_name), collapse = ",")), file = df_snv_class_path, append = T)

# the list of data frames that contain the information of each patient
snv_dfs <- split(df_snv, df_snv$Tumor_Sample_Barcode)

for(i in seq_along(snv_dfs)){
  patient_snv <-left_join(tibble(full_name = unique(df_snv$full_name)), snv_dfs[[i]][!duplicated(snv_dfs[[i]]$full_name), ], by = "full_name")
  # write the information to corresponding csvs
  write(paste0(names(snv_dfs)[i], ",", paste(patient_snv$Variant_Type, collapse = ",")), file = df_snv_type_path, append = T)
  write(paste0(names(snv_dfs)[i], ",", paste(patient_snv$Variant_Classification, collapse = ",")), file = df_snv_class_path, append = T)
}

# delete the downloaded file
unlink(paste0(project_name, "/", "snv_data"), recursive = T)













filename_bla <- list.files(pattern = query_snv[[1]][[1]]$file_name[9], recursive = T)
df_blablabla <- read.maf(filename_bla)
filename_bla <- list.files(pattern = query_snv[[1]][[1]]$file_name[2], recursive = T)
df_blablabla_2 <- read.maf(filename_bla)@data[,c("Hugo_Symbol", "Entrez_Gene_Id", "Gene", 
                                                 "Tumor_Sample_Barcode","Variant_Type", "Variant_Classification")] %>%
  mutate(full_name = paste(Hugo_Symbol, Entrez_Gene_Id, Gene, sep = "|"))













project_name <- "TARGET-AML"
# extracting data from projects start with TARGET-
library(readxl)
# the path to write all the csvs
df_survival_path <- paste0(project_name, "/df_survival.csv")
df_gene_path = paste0(project_name, "/", "df_gene.csv")
df_gene_scale_path = paste0(project_name, "/", "df_gene_scale.csv")
# The path for reading and writing the patient data frames
dir.create(project_name)

# the query to download the survival information of all patients
query_patient_clinical <- try(
  GDCquery(project = project_name,
           data.category = "Clinical",
           data.type = "Clinical Supplement")
)
if(!inherits(query_patient_clinical, "try-error")){
  unlink(paste0(project_name, "/", "patient_data"), recursive = T)
} else {next}
# download the corresponding data
GDCdownload(query_patient_clinical, method = "api", directory = paste0(project_name, "/", "patient_data"))
# extract the filenames that we need
filenames_survival <- query_patient_clinical[[1]][[1]]$file_name[!grepl(pattern = "CDE|Supplement", query_patient_clinical[[1]][[1]]$file_name)]
# the list to store all the survival dfs
dfs_survival <- list()
# loop through all the survival dfs and combine them together
for(i in seq_along(filenames_survival)){
  # acquire the file name of individual .xlsx file
  filename_xlsx <- list.files(pattern = filenames_survival[i], recursive = T)
  # read the .xlsx file
  df_xlsx <- try(
    read_xlsx(filename_xlsx) %>%
    select(patient_id = `TARGET USI`, censoring_status = `Vital Status`, survival_days = `Overall Survival Time in Days`) %>%
    filter(!is.na(censoring_status)) %>%
    mutate(censoring_status = ifelse(censoring_status=="Alive", 0 , 1))
  )
  if(inherits(df_xlsx, "try-error")){next}
  # store them in the list
  dfs_survival[[i]] <- df_xlsx
}
# combine them together
df_survival <- try(rbindlist(dfs_survival))
if(!inherits(df_survival, "try-error")){
  unlink(df_survival_path, recursive = T)
} else {next}
# delete duplicated patient rows
df_survival <- df_survival[!duplicated(df_survival$patient_id), ]

# output the df_survival
fwrite(df_survival, file = df_survival_path)









# project starts with CGCI
library(XML)
data <- xmlParse(file = file.choose())
xml_data <- xmlToList(data)

# projects start with HCMI-CMDC
library(jsonlite)
df_json <- 
  fromJSON(txt = file.choose()) 
# target-os and target-all-p2 need to be fixed