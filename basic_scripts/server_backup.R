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

setwd("~/ShinyApps/project_data")

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

for(i in seq_along(project_id)){
  # the name of the cancer project
  project_name <- project_id[i]  # for example, "TCGA-LUAD"
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
  
  # Create the table for id conversion --------------------------------------
  if(i==1){
    # create individual tables using org.Hs
    egENS <- toTable(org.Hs.egENSEMBL)
    egSYMBOL <- toTable(org.Hs.egSYMBOL)
    
    # bind the tables
    id_table <- egENS %>% left_join(egSYMBOL, by = "gene_id")
    
    # extract the gene ids from df_gene and find their names from id conversion table
    gene_ids <- as_tibble_col(sample_FPKM_name, column_name = "ensembl_id")
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
  # generate the correct name of df_gene and df_gene_scale by joining them with the table that contains all the names
  gene_full_names <-
    left_join(as_tibble_col(gsub("[.]\\d+$","",sample_FPKM_name), column_name = "ensembl_id"),
              full_name_table,
              by = "ensembl_id")$full_name
  # write the gene names as the first line
  write(paste0("patient_id,", paste(gene_full_names, collapse = ",")), file = df_gene_path, append = T)
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
  
  
  # generate the df_gene_scale.csv ------------------------------------------
  # df_gene_scale <- read_csv(df_gene_path)
  df_gene <- data.table::fread(df_gene_path)
  # 
  #   if(i==1){
  #       # create individual tables using org.Hs
  #       egENS <- toTable(org.Hs.egENSEMBL)
  #       egSYMBOL <- toTable(org.Hs.egSYMBOL)
  # 
  #       # bind the tables
  #       id_table <- egENS %>% left_join(egSYMBOL, by = "gene_id")
  # 
  #       # extract the gene ids from df_gene and find their names from id conversion table
  #       gene_ids <- as_tibble_col(colnames(df_gene)[-1], column_name = "ensembl_id")
  #       gene_ids_table <- left_join(gene_ids, id_table, by="ensembl_id")
  # 
  #       # extract the ids that do not exist in the id conversion table
  #       gene_na_table <- gene_ids_table %>%
  #         filter(is.na(gene_id)) %>%
  #         mutate(full_name = ensembl_id)
  #       # extract the ids that exist in the id conversion table
  #       gene_not_na <- gene_ids_table %>%
  #         filter(!is.na(gene_id))
  #       # summarise the duplicated gene_ids and symbols in the df
  #       unique_ids_symbols <-
  #         aggregate(x= tibble(gene_not_na["gene_id"], gene_not_na["symbol"]),
  #                   by=gene_not_na["ensembl_id"],
  #                   FUN = function(X) paste(unique(X), collapse="|")) %>%
  #         mutate(full_name = paste(symbol, gene_id, ensembl_id, sep = "|"))
  #       # the table that contains all the full name of "ENSG...."
  #       full_name_table <- rbind(unique_ids_symbols, gene_na_table)
  #   }
  
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
  # # generate the correct name of df_gene_scale by joining them with the table that contains all the names
  # colnames(df_gene_scale)[-1] <-
  #   left_join(as_tibble_col(colnames(df_gene_scale)[-1], column_name = "ensembl_id"),
  #             full_name_table,
  #             by = "ensembl_id")$full_name
  # delete the previous ones
  if(!inherits(df_gene_scale, "try-error")){
    unlink(df_gene_scale_path, recursive = T)
  } else {next}
  # write the data out
  fwrite(df_gene_scale, file = df_gene_scale_path)
  
  # delete unnecessary files
  unlink(paste0(project_name, "/", "gene_data"), recursive = T)
  unlink(paste0(project_name, "/", "patient_data"), recursive = T)
  
  
  # Create the snv data frames ----------------------------------------------
  
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
  # delete unnecessary files
  unlink(paste0(project_name, "/", "snv_data"), recursive = T)
}


# -------------------------section for TARGET- project------------------------------------------------
# the function to generate df_gene.csv and df_gene_scale.csv
# input: a list of project names; 
# input: shorten_position: the length of the string of case name that you would like to have
# output: two csvs
generate_gene_dfs <- function(project_id, shorten_position = 16L){
  for(i in seq_along(project_id)){
    # the name of the cancer project
    project_name <- project_id[i]  # for example, "TCGA-LUAD"
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
    
    # Create the table for id conversion --------------------------------------
    if(i==1){
      # create individual tables using org.Hs
      egENS <- toTable(org.Hs.egENSEMBL)
      egSYMBOL <- toTable(org.Hs.egSYMBOL)
      
      # bind the tables
      id_table <- egENS %>% left_join(egSYMBOL, by = "gene_id")
      
      # extract the gene ids from df_gene and find their names from id conversion table
      gene_ids <- as_tibble_col(sample_FPKM_name, column_name = "ensembl_id")
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
    # generate the correct name of df_gene and df_gene_scale by joining them with the table that contains all the names
    gene_full_names <-
      left_join(as_tibble_col(gsub("[.]\\d+$","",sample_FPKM_name), column_name = "ensembl_id"),
                full_name_table,
                by = "ensembl_id")$full_name
    # write the gene names as the first line
    write(paste0("patient_id,", paste(gene_full_names, collapse = ",")), file = df_gene_path, append = T)
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
    
    # generate the df_gene_scale.csv ------------------------------------------
    # df_gene_scale <- read_csv(df_gene_path)
    df_gene <- data.table::fread(df_gene_path)
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
  
  # delete unnessary files
  unlink(paste0(project_name, "/", "patient_data"), recursive = T)
  
  
}

#--------generate the .snv files-----------------------------------------
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
  
}

#-------- delete unnecessary files-------------------
unlink("MANIFEST.txt")




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

