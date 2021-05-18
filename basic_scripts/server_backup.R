library(tidyverse)
library(RCurl) # to download the specific folder
library(tools) # to check the file name extensions
library(R.utils) # unzip the .gz files
library(survival) # to do the survival analysis
library(survminer) # to plot the survival analysis nicer
library(data.table)
library(TCGAbiolinks)
# library(microbenchmark) # test the time spent on a command
library(fgsea) # use the function gmtpathways
library("org.Hs.eg.db") # generate the id conversion table
library(maftools) # read the .maf files
library(readxl) # read the .xlsx files
library(TCGAutils) # aliquot UUID to patient barcode conversion

setwd("~/ShinyApps/project_data")
select <- dplyr::select

#-------small functions--------------
return_good_data <- function(df){
  genes <- colnames(df)
  genes_keep <- sapply(genes, function(x){
    data <- df[[x]]
    yn <- table(data[!is.na(data)])
    if(length(yn) == 1){
      return(NULL)
    }else if((length(yn) == 2 & min(yn)>10) | (length(yn) > 2 & sort(yn)[[2]]>10)){
      return(x)
    }else{
      return(NULL)
    }
  })
  genes_keep <- genes_keep[!sapply(genes_keep, is.null)]
  genes_keep <- names(genes_keep)
  if(length(genes_keep)==0){return(NULL)}else{
    df <- data.frame(matrix(0, ncol = length(genes_keep)+1, nrow = 0))
    colnames(df) <- c("patient_id,",genes_keep)
    return(df)
  }
}
# function to remove all na per column
not_all_na <- function(x) any(!is.na(x))
# function to remove all na per column
not_any_na <- function(x) all(!is.na(x))

scaled_df <- function(df_path,sep=",",na="all",na.rm=F){
  odir <- paste0(dirname(df_path),"/tmp/")
  if(!dir.exists(odir)){dir.create(odir)}
  genes <- fread(df_path, sep=sep, nrows = 0)
  genes <- colnames(genes)[-1]
  df <- fread(df_path, sep=sep, select = "patient_id")
  if(na == "all"){
    func_df <- not_all_na
  }else if(na == "any"){
    func_df <- not_any_na
  }
  unlink("df_cnv_scale.csv")
  for(i in seq_along(genes)){
    gene <- genes[i]
    exps <- fread(df_path, sep=sep, select = c("patient_id",gene))
    exps[[gene]] <- as.numeric(exps[[gene]])
    exps[[gene]] <- (exps[[gene]] - mean(exps[[gene]]))/sd(exps[[gene]],na.rm = na.rm)
    if(func_df(exps[[gene]])){
      # # df <- df %>% dplyr::left_join(exps,by="patient_id")
      # fwrite(exps[,2], file = paste0(odir,i,".csv"))
      fwrite(as.list(gene),file="df_cnv_scale.csv", sep = ",", append = T)
    }
  }
  # return(df)
}

#-------get all the project names--------------------
df_gdc_projects <- getGDCprojects()
# write them into a .csv
fwrite(dplyr::select(df_gdc_projects, name, project_id), file = "project_name.csv")
# save the information that are useful to a list
gdc_projects <-list()
for(i in seq_along(df_gdc_projects$project_id)){
  tcga_info <- list()
  tcga_info[[df_gdc_projects$name[i]]][["site"]] <- df_gdc_projects$primary_site[i]
  tcga_info[[df_gdc_projects$name[i]]][["disease"]] <- df_gdc_projects$disease_type[i]
  gdc_projects[[df_gdc_projects$project_id[i]]]<- tcga_info
}

# extract all the tcga projects
df_tcga_projects <- df_gdc_projects %>%
  filter(grepl("TCGA", id))
# extract tcga names to loop through
project_id <- df_tcga_projects$project_id

# function to generate the table containing the full name of all genes using the id conversion table
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
# loop through all tcga projects
# generate the survival dfs
for(w in seq_along(project_id)){
  # the name of the cancer project
  project_name <- project_id[w]  # for example, "TCGA-LUAD"
  # generate the survival data and gene data-------------------------------------
  # the path to write all the csvs
  df_survival_path <- paste0(project_name, "/df_survival.csv")
  # The path for reading and writing the patient data frames
  dir.create(project_name)
  
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
  # create the df_survival
  df_survival <- try(
    read_tsv(filename_patient, skip = 1, na = c("", "NA", "[Not Available]", "[Not Applicable]", "[Discrepancy]")) %>%
      select(patient_id = `bcr_patient_barcode`,gender, `person_neoplasm_cancer_status`, `new_tumor_event_after_initial_treatment`, death_days = `days_to_death`, followup_days = `days_to_last_followup`, diagnosis_days = `days_to_initial_pathologic_diagnosis`) %>%
      mutate(censoring_status = ifelse(is.na(death_days), 0, 1)) %>% # The status indicator, normally 0=alive, 1=dead.
      filter(diagnosis_days == 0) %>%
      mutate(survival_days = as.numeric(ifelse(is.na(death_days), followup_days, death_days))) %>%
      filter(!is.na(survival_days)) %>%
      filter(survival_days >= 0)
  )
  if(!inherits(df_survival, "try-error")){
    unlink(df_survival_path, recursive = T)
  } else {next}
  
  # output the df_survival
  fwrite(df_survival, file = df_survival_path)
  # delete folder after use
  unlink(paste0(project_name, "/", "patient_data"), recursive = T)
  
}


# generate df_gene.csv and df_gene_scale.csv
# loop through all tcga projects
for(w in seq_along(project_id)){
  # the name of the cancer project
  project_name <- project_id[w]  # for example, "TCGA-LUAD"
  df_gene_path = paste0(project_name, "/", "df_gene.csv")
  df_gene_scale_path = paste0(project_name, "/", "df_gene_scale.csv")
  
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
  
  # the name of all cases
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
  sample_FPKM_name <- try(read_tsv(file = gzfile(filename_first_FPKM), col_names = FALSE)$X1)
  if(!inherits(sample_FPKM_name, "try-error")){
    unlink(df_gene_path, recursive = T)
  } else {next}
  
  # write the gene names as the first line
  write(paste0("patient_id,", paste(gsub("[.]\\d+$","",sample_FPKM_name), collapse = ",")), file = df_gene_path, append = T)
  # go through all files and manipulate the duplicates
  for(name in names(cases_info)){
    if(length(cases_info[[name]]) == 1){
      # when no duplicate samples of the patient, write a line to the output file
      filename_gene_info = list.files(pattern = cases_info[[name]][[1]], recursive = T)
      df_gene_write <- read_tsv(file = gzfile(filename_gene_info[1]), col_names = FALSE)
      df_gene_write <- left_join(tibble(X1 = sample_FPKM_name), df_gene_write, by = "X1")
      write(paste0(name, ",", paste(df_gene_write$X2, collapse = ",")), file = df_gene_path, append = T)
      
    } else {
      # the list to store data frames with duplicated patient ids
      df_duplicated_list <- list()
      for(i in seq_along(cases_info[[name]])){
        df_duplicated_list[[i]] <-
          read_tsv(file = gzfile(list.files(pattern = cases_info[[name]][[i]], recursive = T)[1]), col_names = FALSE) %>%
          mutate(X2 = log1p(as.numeric(X2))) # log(x+1)
        names(df_duplicated_list[[i]])[names(df_duplicated_list[[i]]) == "X2"] <- paste0(name, "_", i)
        # print(list.files(pattern = cases_info[[name]][[i]], recursive = T)[1])
      }
      # calculate the means of all duplicates
      df_duplicated <- left_join(tibble(X1 = sample_FPKM_name),
                                 Reduce(function(...) base::merge(..., by = "X1", all.x = T), df_duplicated_list),
                                 by = "X1") %>%
        mutate(mean = rowMeans(select(.,-X1),na.rm = F)) %>%
        mutate(output_col = expm1(mean)) # exp()-1
      write(paste0(name, ",", paste(df_duplicated$output_col, collapse = ",")), file = df_gene_path, append = T)
      
      
    }
  }
  
  
  # generate the df_gene_scale.csv ------------------------------------------
  # df_gene_scale <- read_csv(df_gene_path)
  df_gene <- data.table::fread(df_gene_path)
  
  if(w==1){
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
    # the table that contains all the full name of "ENSG...."
    full_name_table <- rbind(unique_ids_symbols, gene_na_table)
  }
  
  # # start to do the data transformation of df_gene_scale
  patient_ids <- df_gene$patient_id
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
  # generate the correct name of df_gene_scale by joining them with the table that contains all the names
  colnames(df_gene_scale)[-1] <-
    left_join(as_tibble_col(colnames(df_gene_scale)[-1], column_name = "ensembl_id"),
              full_name_table,
              by = "ensembl_id")$full_name
  # delete the previous ones
  if(!inherits(df_gene_scale, "try-error")){
    unlink(df_gene_scale_path, recursive = T)
  } else {next}
  # write the data out
  fwrite(df_gene_scale, file = df_gene_scale_path)
  
  # delete unnecessary files
  unlink(paste0(project_name, "/", "gene_data"), recursive = T)
}

# Create the snv data frames ----------------------------------------------
# loop through all tcga projects
for(w in seq_along(project_id)){
  # the name of the cancer project
  project_name <- project_id[w]  # for example, "TCGA-LUAD"
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
  
  # get all cases(some do not have data) in the project's snv, add to the df later
  project_cases <- getResults(query_snv, cols = c("cases"))
  all_project_cases <- unique(str_sub(colnames(fread(paste0(project_cases[1], "\n"))), start = 1L, end = 12L))
  
  # find the exact file that we want
  special_names <- c("mutect", "varscan", "somaticsniper", "muse")
  for(name_special in special_names){
    # name_special <- "mutect"
    # the path to output the data frames related to snv
    df_snv_type_path = paste0(project_name, "/df_snv_type_", name_special, ".csv")
    df_snv_class_path = paste0(project_name, "/df_snv_class_", name_special, ".csv")
    # read the maf
    filename_snv <- list.files(pattern = grep(name_special, query_snv[[1]][[1]]$file_name,value = T), recursive = T)[1]
    df_maf <- try(read.maf(maf  = filename_snv
                           ,useAll = F
                           ,vc_nonSyn = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins"
                                          ,"Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"
                                          ,"Splice_Site","Translation_Start_Site"
                                          ,"Silent","Splice_Region","Intron","5'UTR","RNA","3'UTR"
                                          ,"5'Flank","3'Flank","IGR")
    ))
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
    write(paste0("patient_id,", paste(unique(df_snv$full_name), collapse = ",")), file = df_snv_type_path, append = T)
    write(paste0("patient_id,", paste(unique(df_snv$full_name), collapse = ",")), file = df_snv_class_path, append = T)
    
    # the list of data frames that contain the infomation of each patient
    snv_dfs <- split(df_snv, df_snv$Tumor_Sample_Barcode)
    
    for(i in seq_along(snv_dfs)){
      # deal with duplicates, delimited by |
      snv_dfs[[i]] <- aggregate(x = select(snv_dfs[[i]], Variant_Type, Variant_Classification),
                                by = list(full_name = snv_dfs[[i]]$full_name),
                                FUN = function(X){paste(unique(X), collapse = "|")}
      )
      # join with all gene names
      patient_snv <-left_join(tibble(full_name = unique(df_snv$full_name)), snv_dfs[[i]][!duplicated(snv_dfs[[i]]$full_name), ], by = "full_name")
      # write the information to corresponding csvs
      write(paste0(names(snv_dfs)[i], ",", paste(patient_snv$Variant_Type, collapse = ",")), file = df_snv_type_path, append = T)
      write(paste0(names(snv_dfs)[i], ",", paste(patient_snv$Variant_Classification, collapse = ",")), file = df_snv_class_path, append = T)
    }
    
    # the patient does not have data from the data we have downloaded
    patient_na <- setdiff(all_project_cases, names(snv_dfs))
    # a vector of NAs
    na_vec <- rep(NA, length(unique(df_snv$full_name)))
    # print their data as NAs in the CSVs
    for(patient_id in patient_na){
      write(paste0(patient_id, ",", paste(na_vec, collapse = ",")), file = df_snv_type_path, append = T)
      write(paste0(patient_id, ",", paste(na_vec, collapse = ",")), file = df_snv_class_path, append = T)
    }
    
  }
  # delete unnecessary files
  unlink(paste0(project_name, "/", "snv_data"), recursive = T)
}
# generate the cnv data frame----------------------------------------------------------
# loop through all tcga projects
for(w in seq_along(project_id)){
  # the name of the cancer project
  project_name <- project_id[w]  # for example, "TCGA-LUAD"
  # the path to output our cnv df
  df_cnv_path <- paste0(project_name, "/df_cnv.csv")
  # the query to download the data related to .cnv file
  query_cnv <- try(
    GDCquery(project = project_name,
             data.category = "Copy Number Variation",
             data.type = "Gene Level Copy Number Scores",
             experimental.strategy = "Genotyping Array",
             platform = "Affymetrix SNP 6.0")
  )
  if(inherits(query_cnv, "try-error")){next}
  # download the corresponding cnv data
  GDCdownload(query_cnv, method = "api", directory = paste0(project_name, "/", "cnv_data"))
  # read the file and clean up the gene ids
  filename_cnv <- list.files(pattern = query_cnv[[1]][[1]]$file_name[1], recursive = T)
  df_cnv <- try(
    fread(filename_cnv) %>%
      mutate(full_name = gsub("[.]\\d+$","",`Gene Symbol`)) %>%
      dplyr::select(full_name, everything()) %>%
      dplyr::select(-`Gene Symbol`))
  if(!inherits(df_cnv, "try-error")){
    unlink(df_cnv_path, recursive = T)
  } else {next}
  # colnames(df_cnv)[colnames(df) == "Gene Symbol"] <- 'full_name'
  # convert the gene ids into the full names of genes
  full_names_cnv <- left_join(as_tibble_col(df_cnv$full_name, column_name = "full_name"),
                              generate_full_name_table(df_cnv$full_name),
                              by = c("full_name" = "ensembl_id"))$full_name.y
  df_cnv <- df_cnv %>%
    mutate(full_name = paste0(full_names_cnv, "|", df_cnv$Cytoband))
  
  # all aliquot UUIDs
  uuids <- colnames(df_cnv)[-3:-1]
  
  # convert to patient barcodes
  patient_ids <- UUIDtoBarcode(uuids, from_type = "aliquot_ids") %>% .[,2]
  
  # column names
  col_names <- c("patient_id",patient_ids)
  
  # remove useless Gene ID and Cytoband columns
  df_cnv <- df_cnv[,c(-2,-3)]
  colnames(df_cnv) <- col_names
  
  
  # deal with duplicate patient ids
  # abbreviate all the patient id into the first 12 characters
  unique_patient_ids <- map_chr(col_names[-1], function(x){
    str_sub(x, start = 1L, end = 12L)
  })
  
  # the list to store all the duplicated and unique patient id information
  duplicated_patient_info <- list()
  for(i in seq_along(unique_patient_ids)){
    patient_shorten_id <- unique_patient_ids[i]
    patient_full_id <- col_names[-1][i]
    if(is.null(duplicated_patient_info[[patient_shorten_id]])){
      duplicated_patient_info[[patient_shorten_id]] <- c(patient_full_id)
    } else {
      duplicated_patient_info[[patient_shorten_id]][length(duplicated_patient_info[[patient_shorten_id]])+1] <- patient_full_id
      # print(patient_shorten_id)
    }
  }
  # write all the gene names as the first line
  write(paste0(colnames(df_cnv)[[1]], ",", paste(df_cnv[[1]], collapse = ",")), file = df_cnv_path, append = T)
  # loop through the list that store duplicated and not duplicated patient ids
  for(patient_shorten_id in names(duplicated_patient_info)){
    if(length(duplicated_patient_info[[patient_shorten_id]]) == 1){
      # if not duplicated, write it out
      write(paste0(patient_shorten_id, ",", paste(df_cnv[[duplicated_patient_info[[patient_shorten_id]][[1]]]], collapse = ",")), file = df_cnv_path, append = T)
    } else {
      # if duplicated, deal with different value
      df_duplicated_patient_info <- df_cnv %>%
        select(duplicated_patient_info[[patient_shorten_id]])
      # the vector that store the resulting values
      output_chr <- c()
      for(n in 1:nrow(df_duplicated_patient_info)){
        if(min(df_duplicated_patient_info[n, ]) >= 0){
          output_chr[n] <- max(df_duplicated_patient_info[n, ])
        } else if(max(df_duplicated_patient_info[n, ]) <= 0){
          output_chr[n] <- min(df_duplicated_patient_info[n, ])
        } else {
          output_chr[n] <- median(as.matrix(df_duplicated_patient_info[n, ]))
        }
      }
      # output it to the csv
      if(!all(output_chr == 0)){
        write(paste0(patient_shorten_id, ",", paste(output_chr, collapse = ",")), file = df_cnv_path, append = T)}
      print(patient_shorten_id)
    }
  }
  # the path to output the df_cnv_scale.csv
  df_cnv_scale_path = paste0(project_name, "/", "df_cnv_scale.csv")
  # read df_cnv into a df
  df_cnv <- try(fread(df_cnv_path))
  # # scale the df
  # df_cnv_scale <- try(
  #   apply(df_cnv[,-1], 2, scale) %>%
  #     as_tibble() %>%
  #     mutate(patient_id = df_cnv$patient_id) %>%
  #     select(patient_id, everything())
  # )
  if(!inherits(df_cnv, "try-error")){
    df_cnv_scale <- return_good_data(df_cnv)
    unlink(df_cnv_scale_path, recursive = T)
    if(!is.null(df_cnv_scale)){
      # output the df to directory
      fwrite(df_cnv_scale, file = df_cnv_scale_path)
    }
  } else {next}
  
  # # delete the cnv data folder
  # unlink(paste0(project_name, "/", "cnv_data"), recursive = T)
}

# --------------------------generate the miRNA data--------------------------------------
# generate the miRNA data
# input: a vector of project id, the position of patient ids that you would like to abbreviate
# output: write to CSVs in the project directory
generate_mir_data <- function(project_ids, abbreviate_position = 12L){
  for(project_name in project_ids){
    # the path to output your data frames
    df_mir_path = paste0(project_name, "/", "df_mir.csv")
    df_mir_scale_path = paste0(project_name, "/", "df_mir_scale.csv")
    
    # the query to download the miRNA information of all patients
    query_mir <- try(
      GDCquery(project = project_name,
               data.category = "Transcriptome Profiling",
               data.type = "miRNA Expression Quantification",
               experimental.strategy = "miRNA-Seq",
               access = "open")
    )
    if((!inherits(query_mir, "try-error")) && !is.null(query_mir)){
      unlink(paste0(project_name, "/", "mir_data"), recursive = T)
    } else {next}
    GDCdownload(query_mir, method = "api", directory = paste0(project_name, "/", "mir_data"))
    # all cases extracted from the query
    cases_name_complete <- query_mir[[1]][[1]]$cases
    # the vector contains the shorten name of all patient ids
    cases_name_shorten <-
      map_chr(cases_name_complete, function(x){
        str_sub(x, start = 1L, end = abbreviate_position)
      })
    # the list to store the duplicated and non-duplicated patient information
    cases_info <- list()
    for(i in seq_along(cases_name_shorten)){
      patient_name <- cases_name_shorten[i]
      file_name <- query_mir[[1]][[1]]$file_name[i]
      if(is.null(cases_info[[patient_name]])){
        # create a new list of unexisted patient
        cases_info[[patient_name]] <- list(file_name)
      } else {
        # store the filename in the list
        cases_info[[patient_name]][[length(cases_info[[patient_name]])+1]] <- file_name
      }
    }
    
    
    # write the gene names into a file a file as the first line
    filename_first_file <- list.files(pattern = cases_info[[1]][[1]], recursive = TRUE)[1]
    sample_mirna_id <- try(fread(filename_first_file)$miRNA_ID)
    if(!inherits(sample_mirna_id, "try-error")){
      unlink(df_mir_path, recursive = T)
    } else {next}
    
    write(paste0("patient_id,", paste(gsub("[.]\\d+$","",sample_mirna_id), collapse = ",")), file = df_mir_path, append = T)
    # go through all files and manipulate the duplicates
    for(name in names(cases_info)){
      if(length(cases_info[[name]]) == 1){
        # when no duplicate samples of the patient, write a line to the output file
        filename_rna_info = list.files(pattern = cases_info[[name]][[1]], recursive = T)
        df_gene_write <- fread(file = filename_rna_info[1]) %>%
          select(miRNA_ID, reads_per_million_miRNA_mapped)
        df_gene_write <- left_join(tibble(miRNA_ID = sample_mirna_id), df_gene_write, by = "miRNA_ID")
        write(paste0(name, ",", paste(df_gene_write$reads_per_million_miRNA_mapped, collapse = ",")), file = df_mir_path, append = T)
        
      } else {
        # the list to store data frames with duplicated patient ids
        df_duplicated_list <- list()
        for(i in seq_along(cases_info[[name]])){
          df_duplicated_list[[i]] <- 
            fread(file = list.files(pattern = cases_info[[name]][[i]], recursive = T)[1]) %>%
            select(miRNA_ID, reads_per_million_miRNA_mapped) %>%
            mutate(reads_per_million_miRNA_mapped = log1p(as.numeric(reads_per_million_miRNA_mapped))) #log(1+x)
          names(df_duplicated_list[[i]])[names(df_duplicated_list[[i]]) == "reads_per_million_miRNA_mapped"] <- paste0(name, "_", i)
          # print(list.files(pattern = cases_info[[name]][[i]], recursive = T)[1])
        }
        # calculate the means of all duplicates
        df_duplicated <- left_join(tibble(miRNA_ID = sample_mirna_id),
                                   Reduce(function(...) base::merge(..., by = "miRNA_ID", all.x = T), df_duplicated_list),
                                   by = "miRNA_ID") %>%
          mutate(mean = rowMeans(select(.,-miRNA_ID),na.rm = F)) %>%
          mutate(output_col = expm1(mean)) # exp()-1
        write(paste0(name, ",", paste(df_duplicated$output_col, collapse = ",")), file = df_mir_path, append = T)
        
        
      }
    }
    
    # read df_mir into a df
    df_mir <- fread(df_mir_path)
    # scale the df
    df_mir_scale <- try(
      apply(df_mir[,-1], 2, scale) %>%
        as_tibble() %>%
        mutate(patient_id = df_mir$patient_id) %>%
        select(patient_id, everything())
    )
    if(!inherits(df_mir_scale, "try-error")){
      unlink(df_mir_scale_path, recursive = T)
    }
    # output the df to directory
    fwrite(df_mir_scale, file = df_mir_scale_path)
    unlink(paste0(project_name, "/", "mir_data"), recursive = T)
  }
}

# apply the function to generate miRNA data for TCGA projects
generate_mir_data(project_ids = project_id)

# generate the methylation data----------------------------------------------------

# the function to get all the gene names in the format of "cg00013618|IGLVI-70|PRAMENP"
# input: a df tthat from met data
# output: a vector of all unique gene names
get_met_names <- function(df){
  unique_gene_symbols <- c()
  for(i in seq_along(df$Gene_Symbol)){
    x <- df$Gene_Symbol[i]
    cg_id <- df$`Composite Element REF`[i]
    symbol <- paste(unique(base::strsplit(x, split = ";")[[1]]), collapse = "|")
    if(symbol == "."){unique_gene_symbols[i] <- cg_id}
    else{unique_gene_symbols[i] <- paste(cg_id, symbol, sep = "|")}
  }
  unique_gene_symbols
}

# the function to scale and write a data frame, for example, df_mir
# input: a path of original df; a path to output the scaled df
# output the df to the path you specify as the second argument
scale_a_df <- function(df_mir_path, df_mir_scale_path){
  # read df_mir into a df
  df_mir <- fread(df_mir_path)
  # scale the df
  df_mir_scale <- try(
    apply(df_mir[,-1], 2, scale) %>%
      as_tibble() %>%
      mutate(patient_id = df_mir$patient_id) %>%
      select(patient_id, everything())
  )
  if(!inherits(df_mir_scale, "try-error")){
    unlink(df_mir_scale_path, recursive = T)
  }
  # output the df to directory
  fwrite(df_mir_scale, file = df_mir_scale_path)
}

# start generating csvs-------------------------------------------------------
for(w in seq_along(project_id)){
  # the name of the cancer project
  project_name <- project_id[w]  # for example, "TCGA-LUAD"
  # the path to output your data frames
  df_met_path = paste0(project_name, "/", "df_met.csv")
  df_met_scale_path = paste0(project_name, "/", "df_met_scale.csv")
  # the query to download the methylation information of all patients
  query_met <- try(
    GDCquery(project = project_name,
             data.category = "DNA Methylation",
             data.type = "Methylation Beta Value",
             experimental.strategy = "Methylation Array",
             platform = "Illumina Human Methylation 450")
  )
  if(!inherits(query_met, "try-error")){
    unlink(paste0(project_name, "/", "met_data"), recursive = T)
  } else {next}
  GDCdownload(query_met, method = "api", directory = paste0(project_name, "/", "met_data"))
  # all cases extracted from the query
  cases_name_complete <- query_met[[1]][[1]]$cases
  # the vector contains the shorten name of all patient ids
  cases_name_shorten <-
    map_chr(cases_name_complete, function(x){
      str_sub(x, start = 1L, end = 12L)
    })
  # the list to store the duplicated and non-duplicated patient information
  cases_info <- list()
  for(i in seq_along(cases_name_shorten)){
    patient_name <- cases_name_shorten[i]
    file_name <- query_met[[1]][[1]]$file_name[i]
    if(is.null(cases_info[[patient_name]])){
      # create a new list of unexisted patient
      cases_info[[patient_name]] <- list(file_name)
    } else {
      # store the filename in the list
      cases_info[[patient_name]][[length(cases_info[[patient_name]])+1]] <- file_name
    }
  }
  
  
  # write the gene names into a file a file as the first line
  filename_first_file <- list.files(pattern = cases_info[[1]][[1]], recursive = TRUE)[1]
  unique_gene_names <- try(get_met_names(df = fread(file = filename_first_file)))
  if(!inherits(unique_gene_names, "try-error")){
    unlink(df_met_path, recursive = T)
  } else {next}
  
  write(paste0("patient_id,", paste(unique_gene_names, collapse = ",")), file = df_met_path, append = T)
  # go through all files and manipulate the duplicates
  for(name in names(cases_info)){
    if(length(cases_info[[name]]) == 1){
      # when no duplicate samples of the patient, write a line to the output file
      filename_met_info = list.files(pattern = cases_info[[name]][[1]], recursive = T)
      write(paste0(name, ",", paste(fread(file = filename_met_info[1])$Beta_value, collapse = ",")), file = df_met_path, append = T)
      
    } else {
      # the list to store data frames with duplicated patient ids
      df_duplicated_list <- list()
      for(i in seq_along(cases_info[[name]])){
        df_duplicated_list[[i]] <- 
          fread(file = list.files(pattern = cases_info[[name]][[i]], recursive = T)[1]) %>%
          select(`Composite Element REF`, Beta_value) %>%
          mutate(Beta_value = log1p(as.numeric(Beta_value))) #log(1+x)
        # assign name to avoid duplicated column names
        names(df_duplicated_list[[i]])[names(df_duplicated_list[[i]]) == "Beta_value"] <- paste0(name, "_", i)
        # print(list.files(pattern = cases_info[[name]][[i]], recursive = T)[1])
      }
      # calculate the means of all duplicates
      df_duplicated <- Reduce(function(...) base::merge(..., by = "Composite Element REF", all.x = T), df_duplicated_list) %>%
        mutate(mean = rowMeans(select(., -"Composite Element REF"),na.rm = F)) %>%
        mutate(output_col = expm1(mean)) # exp()-1
      write(paste0(name, ",", paste(df_duplicated$output_col, collapse = ",")), file = df_met_path, append = T)
      
      
    }
  }
  
  # scale df_met to df_met_scale
  scale_a_df(df_met_path, df_met_scale_path)
  # delete the folder to save memory
  unlink(x = paste0(project_name, "/", "met_data"), recursive = T)
}

# -------------------------section for TARGET- project------------------------------------------------
# the function to generate df_gene.csv and df_gene_scale.csv
# input: a list of project names; 
# input: shorten_position: the length of the string of case name that you would like to have
# output: two csvs
generate_gene_dfs <- function(project_id, shorten_position = 16L){
  for(k in seq_along(project_id)){
    # the name of the cancer project
    project_name <- project_id[k]  # for example, "TCGA-LUAD"
    # the path to write all the csvs
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
    
    # the vector contain all patient ids
    cases_name_complete <- query_gene_counts[[1]][[1]]$cases
    # the vector contains the shorten name of all patient ids
    cases_name_shorten <-
      map_chr(cases_name_complete, function(x){
        str_sub(x, start = 1L, end = shorten_position)
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
    
    # # Create the table for id conversion --------------------------------------
    # if(i==1){
    #   # create individual tables using org.Hs
    #   egENS <- toTable(org.Hs.egENSEMBL)
    #   egSYMBOL <- toTable(org.Hs.egSYMBOL)
    #   
    #   # bind the tables
    #   id_table <- egENS %>% left_join(egSYMBOL, by = "gene_id")
    #   
    #   # extract the gene ids from df_gene and find their names from id conversion table
    #   gene_ids <- as_tibble_col(sample_FPKM_name, column_name = "ensembl_id")
    #   gene_ids_table <- left_join(gene_ids, id_table, by="ensembl_id")
    #   
    #   # extract the ids that do not exist in the id conversion table
    #   gene_na_table <- gene_ids_table %>%
    #     filter(is.na(gene_id)) %>%
    #     mutate(full_name = ensembl_id)
    #   # extract the ids that exist in the id conversion table
    #   gene_not_na <- gene_ids_table %>%
    #     filter(!is.na(gene_id))
    #   # summarise the duplicated gene_ids and symbols in the df
    #   unique_ids_symbols <-
    #     aggregate(x= tibble(gene_not_na["gene_id"], gene_not_na["symbol"]),
    #               by=gene_not_na["ensembl_id"],
    #               FUN = function(X) paste(unique(X), collapse="|")) %>%
    #     mutate(full_name = paste(symbol, gene_id, ensembl_id, sep = "|"))
    #   # the table that contains all the full name of "ENSG...."
    #   full_name_table <- rbind(unique_ids_symbols, gene_na_table)
    # }
    # # generate the correct name of df_gene and df_gene_scale by joining them with the table that contains all the names
    # gene_full_names <-
    #   left_join(as_tibble_col(gsub("[.]\\d+$","",sample_FPKM_name), column_name = "ensembl_id"),
    #             full_name_table,
    #             by = "ensembl_id")$full_name
    # write the gene names as the first line
    write(paste0("patient_id,", paste(gsub("[.]\\d+$","",sample_FPKM_name), collapse = ",")), file = df_gene_path, append = T)
    # go through all files and manipulate the duplicates
    for(name in names(cases_info)){
      if(length(cases_info[[name]]) == 1){
        # when no duplicate samples of the patient, write a line to the output file
        filename_gene_info = list.files(pattern = cases_info[[name]][[1]], recursive = T)
        df_gene_write <- read_tsv(file = gzfile(filename_gene_info[1]), col_names = FALSE)
        df_gene_write <- left_join(tibble(X1 = sample_FPKM_name), df_gene_write, by = "X1")
        write(paste0(name, ",", paste(df_gene_write$X2, collapse = ",")), file = df_gene_path, append = T)
        
      } else {
        # the list to store data frames with duplicated patient ids
        df_duplicated_list <- list()
        for(i in seq_along(cases_info[[name]])){
          df_duplicated_list[[i]] <-
            read_tsv(file = gzfile(list.files(pattern = cases_info[[name]][[i]], recursive = T)[1]), col_names = FALSE) %>%
            mutate(X2 = log1p(as.numeric(X2))) # log(x+1)
          names(df_duplicated_list[[i]])[names(df_duplicated_list[[i]]) == "X2"] <- paste0(name, "_", i)
          # print(list.files(pattern = cases_info[[name]][[i]], recursive = T)[1])
        }
        # calculate the means of all duplicates
        df_duplicated <- left_join(tibble(X1 = sample_FPKM_name),
                                   Reduce(function(...) base::merge(..., by = "X1", all.x = T), df_duplicated_list),
                                   by = "X1") %>%
          mutate(mean = rowMeans(select(.,-X1),na.rm = F)) %>%
          mutate(output_col = expm1(mean)) # exp()-1
        write(paste0(name, ",", paste(df_duplicated$output_col, collapse = ",")), file = df_gene_path, append = T)
        
        
      }
    }
    
    # generate the df_gene_scale.csv ------------------------------------------
    # df_gene_scale <- read_csv(df_gene_path)
    df_gene <- data.table::fread(df_gene_path)
    # create the table for gene id conversion to full names
    if(k==1){
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
      # the table that contains all the full name of "ENSG...."
      full_name_table <- rbind(unique_ids_symbols, gene_na_table)
    }
    
    # start to do the data transformation of df_gene_scale
    patient_ids <- df_gene$patient_id
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
    # generate the correct name of df_gene_scale by joining them with the table that contains all the names
    colnames(df_gene_scale)[-1] <-
      left_join(as_tibble_col(colnames(df_gene_scale)[-1], column_name = "ensembl_id"),
                full_name_table,
                by = "ensembl_id")$full_name
    
    if(!inherits(df_gene_scale, "try-error")){
      unlink(df_gene_scale_path, recursive = T)
    } else {next}
    # write the data out
    fwrite(df_gene_scale, file = df_gene_scale_path)
    
    # delete unnecessary files
    unlink(paste0(project_name, "/", "gene_data"), recursive = T)
  }
}

# Start generating files for TARGET- files --------------------------------
# acquire all the TARGET projects
df_target_projects <- df_gdc_projects %>%
  filter(grepl("TARGET", id))
project_id <- df_target_projects$project_id
# generate the df_gene.csv and df_gene_scale.csv using the function--------
generate_gene_dfs(project_id)
# loop through all projects
for(j in seq_along(project_id)){
  # ----------------generate the df_survival.csv--------------------------------------
  project_name <- project_id[j]
  # the path to write all the csvs
  df_survival_path <- paste0(project_name, "/df_survival.csv")
  
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
        dplyr::select(patient_id = `TARGET USI`,gender = Gender, censoring_status = `Vital Status`, survival_days = `Overall Survival Time in Days`) %>%
        filter(!is.na(censoring_status)) %>%
        mutate(censoring_status = ifelse(censoring_status=="Alive", 0 , 1)) %>%
        filter(survival_days >= 0)
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
  
  # delete unnessary files
  unlink(paste0(project_name, "/", "patient_data"), recursive = T)
  
  
}

#--------generate the snv files-----------------------------------------
for(k in seq_along(project_id)){
  project_name <- project_id[k]
  print(project_name)
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
  
  # check if the df is empty; if yes, skip to the next iteration
  if(nrow(df_snv) == 0 && ncol(df_snv) == 0){
    unlink(paste0(project_name, "/", "snv_data"), recursive = T)
    unlink(list.files(pattern = ".tar.gz"))
    next
  }
  
  # write the first line of the csvs
  write(paste0("patient_id,", paste(unique(df_snv$full_name), collapse = ",")), file = df_snv_type_path, append = T)
  write(paste0("patient_id,", paste(unique(df_snv$full_name), collapse = ",")), file = df_snv_class_path, append = T)
  
  # the list of data frames that contain the information of each patient
  snv_dfs <- split(df_snv, df_snv$Tumor_Sample_Barcode)
  
  for(i in seq_along(snv_dfs)){
    # deal with duplicates, delimited by |
    snv_dfs[[i]] <- aggregate(x = select(snv_dfs[[i]], Variant_Type, Variant_Classification),
                              by = list(full_name = snv_dfs[[i]]$full_name),
                              FUN = function(X){paste(unique(X), collapse = "|")}
    )
    # join with all gene names
    patient_snv <-left_join(tibble(full_name = unique(df_snv$full_name)), snv_dfs[[i]][!duplicated(snv_dfs[[i]]$full_name), ], by = "full_name")
    # write the information to corresponding csvs
    write(paste0(names(snv_dfs)[i], ",", paste(patient_snv$Variant_Type, collapse = ",")), file = df_snv_type_path, append = T)
    write(paste0(names(snv_dfs)[i], ",", paste(patient_snv$Variant_Classification, collapse = ",")), file = df_snv_class_path, append = T)
  }
  
  # delete the downloaded file
  unlink(paste0(project_name, "/", "snv_data"), recursive = T)
}



# generate the cnv file----------------------------------------
# loop through all projects
for(j in seq_along(project_id)){
  # ----------------generate the df_survival.csv--------------------------------------
  project_name <- project_id[j]
  # the path to output our cnv df
  df_cnv_path <- paste0(project_name, "/df_cnv.csv")
  # the query to download the data related to .cnv file
  query_cnv <- try(
    GDCquery(project = project_name,
             data.category = "Copy Number Variation",
             data.type = "Allele-specific Copy Number Segment")
  )
  if(inherits(query_cnv, "try-error")){next}
  # download the corresponding cnv data
  GDCdownload(query_cnv, method = "api", directory = paste0(project_name, "/", "cnv_data"))
  
  # the abbreviated patient ids
  patient_id_TARGET <- try(
    map_chr(query_cnv[[1]][[1]]$cases, function(x){
      str_sub(x, start = 1L, end = 16L)
    })
  )
  if(!inherits(patient_id_TARGET, "try-error")){
    unlink(df_cnv_path, recursive = T)
  } else {next}
  
  chromosome <- list()
  df_target_cnv <- list()
  # store all dfs and their chromosomes in lists
  for(i in seq_along(query_cnv[[1]][[1]]$file_name)){
    # the df storing information for a patient
    df_single_patient <- fread(list.files(pattern = query_cnv[[1]][[1]]$file_name[i], recursive = T)[1]) %>%
      mutate(full_name = paste0(Chromosome, ":", Start, "-", End)) %>%
      select(full_name, Copy_Number)
    if(is.null(df_target_cnv[[patient_id_TARGET[i]]])){
      df_target_cnv[[patient_id_TARGET[i]]][[1]] <- df_single_patient
    } else{
      df_target_cnv[[patient_id_TARGET[i]]][[length(df_target_cnv[[patient_id_TARGET[i]]]) + 1]] <- df_single_patient
    }
    chromosome[[i]] <- df_single_patient$full_name
  }
  
  # acquire all unique chromosomes in lists
  unique_chromosome <- unique(unlist(chromosome))
  # write the columns names as the first line of the .csv
  write(paste0("patient_id,", paste(unique_chromosome, collapse = ",")), file = df_cnv_path, append = T)
  
  # loop through each unique patient
  for(j in seq_along(df_target_cnv)){
    # the patient has only one case
    if(length(df_target_cnv[[j]]) == 1){
      df_output<- left_join(as_tibble_col(unique_chromosome, column_name = "full_name"), df_target_cnv[[j]][[1]], by = "full_name")
      write(paste0(names(df_target_cnv)[j], ",", paste(df_output$Copy_Number, collapse = ",")), file = df_cnv_path, append = T)
    } else {
      # the patient that has more than one cases, join all data frames
      df_output <- rbindlist(df_target_cnv[[j]], use.names = TRUE)
      # deal with duplicated values
      df_output <- left_join(as_tibble_col(unique_chromosome, column_name = "full_name"),
                             aggregate(x= df_output$Copy_Number,
                                       by=list(full_name = df_output$full_name),
                                       FUN = function(X){
                                         if(max(X) <= 2){
                                           return(min(X))
                                         }
                                         else if(min(X) >=2 ){
                                           return(max(X))
                                         }
                                         else {return(median(X))}
                                       }), by = "full_name")
      # output the results
      if(!all(df_output$x == 2)){
        write(paste0(names(df_target_cnv)[j], ",", paste(df_output$x, collapse = ",")), file = df_cnv_path, append = T)}
    }
  }
  
  # the path to output the df_cnv_scale.csv
  df_cnv_scale_path = paste0(project_name, "/", "df_cnv_scale.csv")
  # read df_cnv into a df
  df_cnv <- try(fread(df_cnv_path))
  # # scale the df
  # df_cnv_scale <- try(
  #   apply(df_cnv[,-1], 2, scale) %>%
  #     as_tibble() %>%
  #     mutate(patient_id = df_cnv$patient_id) %>%
  #     select(patient_id, everything())
  # )
  if(!inherits(df_cnv, "try-error")){
    df_cnv_scale <- return_good_data(df_cnv)
    unlink(df_cnv_scale_path, recursive = T)
    if(!is.null(df_cnv_scale)){
      # output the df to directory
      fwrite(df_cnv_scale, file = df_cnv_scale_path)
    }
  } else {next}
  
  # # delete unneccessary folder
  # unlink(paste0(project_name, "/", "cnv_data"), recursive = T)
  # 
}

# apply the function to generate miRNA data for TARGET projects
generate_mir_data(project_ids = project_id, abbreviate_position = 16L)

#-------- delete unnecessary files-------------------
unlink("MANIFEST.txt")

# generate a data frame that contains all the  data name of TARGET projects
write("project_name,existing_data", file = "TARGET_existing_data.csv")
for(project_name in df_target_projects$project_id){
  existing_data <- c()
  if(length(list.files(path = project_name, pattern = "df_gene_scale")) > 0){
    existing_data[length(existing_data)+1] <- "rna"
  }
  if(length(list.files(path = project_name, pattern = "df_snv")) > 0){
    existing_data[length(existing_data)+1] <- "snv"
  }
  if(length(list.files(path = project_name, pattern = "df_cnv_scale")) > 0){
    existing_data[length(existing_data)+1] <- "cnv"
  }
  if(length(list.files(path = project_name, pattern = "df_mir_scale")) > 0){
    existing_data[length(existing_data)+1] <- "mir"
  }
  write(paste0("\"", project_name, "\",\"", paste(existing_data, collapse = ","), "\"") , file = "TARGET_existing_data.csv", append = T)
}

# generate a data frame that contains all the missing data name of TARGET projects
write("project_name,missing_data", file = "TARGET_missing_data.csv")
for(project_name in df_target_projects$project_id){
  missing_data <- c()
  if(length(list.files(path = project_name, pattern = "df_gene_scale")) == 0){
    missing_data[length(missing_data)+1] <- "rna"
  }
  if(length(list.files(path = project_name, pattern = "df_snv")) == 0){
    missing_data[length(missing_data)+1] <- "snv"
  }
  if(length(list.files(path = project_name, pattern = "df_cnv_scale")) == 0){
    missing_data[length(missing_data)+1] <- "cnv"
  }
  if(length(list.files(path = project_name, pattern = "df_mir_scale")) == 0){
    missing_data[length(missing_data)+1] <- "mir"
  }
  write(paste0("\"", project_name, "\",\"", paste(missing_data, collapse = ","), "\"") , file = "TARGET_missing_data.csv", append = T)
}

# generate a data frame that contains all the missing data names of TCGA projects
write("project_name,missing_data", file = "TCGA_missing_data.csv")
for(project_name in df_tcga_projects$project_id){
  missing_data <- c()
  if(length(list.files(path = project_name, pattern = "df_gene")) == 0){
    missing_data[length(missing_data)+1] <- "rna"
  }
  if(length(list.files(path = project_name, pattern = "df_snv")) == 0){
    missing_data[length(missing_data)+1] <- "snv"
  }
  if(length(list.files(path = project_name, pattern = "df_cnv")) == 0){
    missing_data[length(missing_data)+1] <- "cnv"
  }
  if(length(list.files(path = project_name, pattern = "df_mir")) == 0){
    missing_data[length(missing_data)+1] <- "mir"
  }
  if(length(list.files(path = project_name, pattern = "df_met")) == 0){
    missing_data[length(missing_data)+1] <- "mir"
  }
  write(paste0("\"", project_name, "\",\"", paste(missing_data, collapse = ","), "\"") , file = "TCGA_missing_data.csv", append = T)
}


#--------- test_code-------------------------------#
project_name <- "TARGET-AML"
dir.create(project_name)
# create a query that downloads all the snv data
query_snv <- try(
  GDCquery(project = project_name,
           data.category = "Simple Nucleotide Variation",
           data.type = "Masked Somatic Mutation",
           experimental.strategy = "WXS")
)
if(inherits(query_snv, "try-error")){next}
GDCdownload(query_snv, method = "api", directory = paste0(project_name, "/", "snv_data"))
df_target_snv_list <- list()
for(j in seq_along(query_snv[[1]][[1]]$file_name)){
  filename_target_snv <- list.files(pattern = query_snv[[1]][[1]]$file_name[j], recursive = T)
  df_maf_target_snv <- read.maf(filename_target_snv)@data[,c("Hugo_Symbol", "Entrez_Gene_Id", "Gene", 
                                                             "Tumor_Sample_Barcode","Variant_Type", "Variant_Classification")] %>%
    mutate(full_name = paste(Hugo_Symbol, Entrez_Gene_Id, Gene, sep = "|"))
  df_target_snv_list[[j]] <- df_maf_target_snv
}

