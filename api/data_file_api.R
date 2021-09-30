library(stringr)
library(hash)
library(plumber)
library(tidyverse)
library(data.table)

# source("~/ShinyApps/global/variables.R")

#* @apiTitle cSurvival API
#* @apiDescription API to download cSurvival's database
#* @apiVersion 1.0.0

#* @apiTitle cSurvival API
#* @apiDescription API to download cSurvival's database
#* @apiVersion 1.0.0

api_path = "/home/cSurvival/ShinyApps"
indir <- paste0(api_path, "/www/project_data/")
df_ccle <- fread(paste0(indir,"DepMap/ccle.csv"), select = c("patient_id","CCLE_Name","gender","primary_or_metastasis","primary_disease","Subtype"))

#Initialize a dictionary to record all possible get data combination----
DepMap_dictionary <- hash()
DepMap_dictionary[["Mutation"]] <- "df_snv_class.csv"
DepMap_dictionary[["CNV"]] <- "df_cnv.csv"
DepMap_dictionary[["Proteomics"]] <- "df_proteomic.csv"
DepMap_dictionary[["CRISPR-Cas9"]] <- "DepMap-CRISPR.csv"
DepMap_dictionary[["Expression"]] <- "df_gene.csv"
DepMap_dictionary[["RNAi"]] <- "DepMap-RNAi.csv"
DepMap_dictionary[["Drug"]] <- "DepMap-Drug.csv"

TARGET_dictionary <- hash()
TARGET_dictionary[["Expression"]] <- "df_gene.csv"
TARGET_dictionary[["Clinical"]] <- "df_survival.csv"
TARGET_dictionary[["miRNA"]] <- "df_mir.csv"
TARGET_dictionary[["CNV"]] <- "df_cnv.csv"
TARGET_dictionary[["Mutation"]] <- "df_snv_class.csv"


TCGA_dictionary <- hash()
TCGA_dictionary[["Expression"]] <- "df_gene.csv"
TCGA_dictionary[["Mutation"]] <- "df_snv_class_977.csv"
TCGA_dictionary[["CNV"]] <- "df_cnv.csv"
TCGA_dictionary[["Clinical"]] <- "df_survival_o.csv"
TCGA_dictionary[["miRNA"]] <- "df_mir.csv"
TCGA_dictionary[["Methylation"]] <- "df_met.csv"
TCGA_dictionary[["RPPA"]] <- "df_rrpa.csv"

#helper function for error handler
error_response <- function(res, project_name, error_type){
  message = hash() 
  message[["file_command"]] <- paste0("The file command you entered cannot be found in ",project_name ," API, please check API description again for acceptable file command")
  message[["subproject"]] <- paste0("The subproject you entered cannot be found in  ",project_name," API, please check cSurvival website again for available ",project_name," projects")
  message[["no_such_file"]] <- paste0("This subproject does not have the file you requested, please check cSurvival website again for available files of this subproject")
  message[["no_queried_cancer_types_found"]] <- paste0("No queried cancer types can be found in ", project_name, ".")
  message[["no_queried_cancer_subtypes_found"]] <- paste0("No queried cancer subtypes can be found in ", project_name, ".")
  message[["no_queried_cell_lines_found"]] <- paste0("No queried cell lines can be found in ", project_name, ".")
  message[["queried_gene_not_found"]] <- paste0("The queried gened cannot be found in ", project_name, ".")
  message[["invalid_no_of_analysis"]] <- paste0("Invalid no of analysis. Only no_of_analysis equal to 1 or 2 is supported.")
  res$status <- 500
  res$body <- jsonlite::toJSON(auto_unbox = TRUE, list(
    status = 500,
    message = message[[error_type]]
  ))
  return(res)
}



#* Get specific data inside of DepMap projects, acceptable file names are: "Expression", "Mutation","CNV" ,"Proteomics", "CRISPR-Cas9", "RNAi", "Drug"
#* @param file:[str] subproject name_file name
#* @serializer contentType list(type="text/csv")
#* @get /DepMap_Data/<file>
function(file = "Expression", res){
  # #Set work directory to current file
  # setwd(api_path)
  #!!! subject to change
  #handle file command error
  if(!(file %in% keys(DepMap_dictionary))){
    return(error_response(res, project_name = "DepMap", error_type = "file_command"))
  }
  #and filename(to get the exact file based on query command)
  filename = DepMap_dictionary[[file]]
  filename <- paste0(api_path,"/www/project_data/DepMap/",filename)
  cat("The directory api is trying to access is: ",filename,"\n")
  
  if(!(file.exists(filename))){
    return(error_response(res, project_name = "DepMap", error_type = "no_such_file"))
  }
  plumber::include_file(filename, res, "text/csv")
}



#* Get specific data inside of TCGA projects, acceptable file names are: "Expression", "Mutation","CNV", "Clinical", "miRNA", "Methylation", "RPPA"
#* @param file:[str] subproject name_file name
#* @param file:[str] subproject name  + file name
#* @serializer contentType list(type="text/csv")
#* @get /TCGA_Data/<file>
function(file = "LUSC_Mutation", res){
  #Set work directory to current file
  setwd("/Applications/Codes/cSurvival_dev")
  #!!! subject to change
  #split the query into two parts: subproject(folder)
  subproject = unlist(str_split(file, "_", n = 2))[1]
  #and filename(to get the exact file based on query command)
  filename = unlist(str_split(file, "_", n = 2))[2]
  #handle file command error
  if(!(filename %in% keys(TCGA_dictionary))){
    return(error_response(res, project_name = "TCGA", error_type = "file_command"))
  }
  #handle subproject error
  if(!(file.exists(paste0(api_path, "/www/project_data/TCGA-",subproject)))){
    return(error_response(res, project_name = "TCGA", error_type = "subproject"))
  }
  
  filename = TCGA_dictionary[[filename]]
  filename <- paste0(api_path, "/www/project_data/","TCGA-",subproject,"/",filename)
  cat("The directory api is trying to access is: ",filename,"\n")
  if(!(file.exists(filename))){
    return(error_response(res, project_name = "TCGA", error_type = "no_such_file"))
  }
  plumber::include_file(filename, res, "text/csv")
}


#* Get specific data inside of TARGET projects, acceptable file commands are: "Expression", "Mutation" , "Clinical", "miRNA", "CNV"
#* @param file:[str] subproject name_file command
#####* @serializer contentType list(type="text/csv")
#* @get /TARGET_Data/<file>
function(file = "ALL-P1_Mutation", res){
  #Set work directory to current file
  setwd("/Applications/Codes/cSurvival_dev")
  #!!! subject to change
  #split the query into two parts: subproject(folder)
  subproject = unlist(str_split(file, "_", n = 2))[1]
  #and filename(to get the exact file based on query command)
  filename = unlist(str_split(file, "_", n = 2))[2]
  
  #handle file command error
  if(!(filename %in% keys(TARGET_dictionary))){
    return(error_response(res, project_name = "TARGET", error_type = "file_command"))
  }
  #handle subproject error
  if(!(file.exists(paste0(api_path, "/www/project_data/TARGET-",subproject)))){
    return(error_response(res, project_name = "TARGET", error_type = "subproject"))
  }
  
  
  filename = TARGET_dictionary[[filename]]
  filename <- paste0(api_path, "/www/project_data/","TARGET-",subproject,"/",filename)
  cat("The directory api is trying to access is: ",filename,"\n")
  #handle specific file not found
  if(!(file.exists(filename))){
    return(error_response(res, project_name = "TARGET", error_type = "no_such_file"))
  }
  plumber::include_file(filename, res, "text/csv")
}


#* Some descrption here
#* @param analysis_parameter e.g project=&no_of_analysis=&cancer_type=&cancer_subtype=&cell_line=&gene1=&data_category&=&data_type=&gene2=
#* @get /DepMap_Analysis
function(analysis_parameter, res){
  parameter_list <- unlist(str_split(analysis_parameter,"&"))
  parameter_dict <- hash()
  subtypes_optional_default_all <- T
  cell_lines_optional_default_all <- T

  for (par in parameter_list) {
    curr <- unlist(str_split(par,"="))
    parameter_dict[[curr[1]]] <- curr[2]
  }

  # user inputs read from api command
  project <- parameter_dict[["project"]]
  variable_n <- parameter_dict[["no_of_analysis"]]
  # multiple cancer types/subtypes/cell lines seperated by semicolon are allowed
  cancer_type <- parameter_dict[["cancer_type"]]
  ccle_cancer_types <- unlist(str_split(cancer_type,";"))
  if ("cancer_subtype" %in% keys(parameter_dict)) {
    cancer_subtype <- parameter_dict[["cancer_subtype"]]
    ccle_cancer_subtypes <- unlist(str_split(cancer_subtype,";"))
    subtypes_optional_default_all <- F
  }
  if ("cell_line" %in% keys(parameter_dict)) {
    cell_line <- parameter_dict[["cell_line"]]
    cell_line <- unlist(str_split(cell_line,";"))
    cell_lines_optional_default_all <- F
  }

  # if (variable_n > 2 | variable_n < 1) {
  #   return (error_response(res, project_name = "DepMap", error_type="invalid_no_of_analysis"))
  # } else {
  #   if (variable_n == 1) {
  #     gene1 <- parameter_dict[["gene1"]]
  #     data_category1 <- parameter_dict[["data_category1"]]
  #   }
  # }

  plot_stype <- paste0(gsub("^DepMap-","",project)," dependency")
  depmap_path <- paste0(indir,"DepMap/DepMap-",project,".csv")
  # the genes for initial selection
  depmap_genes <- fread(depmap_path, sep = ",", nrows = 0, quote="")
  depmap_genes <- colnames(depmap_genes)[-1]
  # available cell lines
  depmap_ids <- fread(depmap_path, sep = ",", select = "patient_id", quote="") %>% .[["patient_id"]]
  depmap_ccle <- df_ccle %>% dplyr::filter(patient_id %in% depmap_ids)
  cell_lines <- depmap_ccle$CCLE_Name
  names(cell_lines) <- depmap_ccle$patient_id

  # primary cancers - options available in "Select cancer type(s) to study"
  ccle_cancers <- sort(unique(depmap_ccle$primary_disease))
  # check whether cancer types entered are valid or not
  cancer_not_found <- list()
  cancer_found <- list()
  for (cancer in ccle_cancer_types) {
    cancer <- str_trim(cancer)
    if (tolower(cancer) %in% tolower(ccle_cancers)) {
      cancer_found <- append(cancer_found, cancer)
    } else {
      cancer_not_found <- append(cancer_not_found, cancer)
    }
  }
  if (length(cancer_found) == 0) {
    return (error_response(res, project_name = "DepMap", error_type="no_queried_cancer_types_found"))
  } else {
    ccle_cancer_types <- cancer_found
    # process cancer types selected to get all the corresponding subtypes
    subtypes <- depmap_ccle %>% dplyr::filter(tolower(primary_disease) %in% tolower(ccle_cancer_types)) %>%
      .[["Subtype"]] %>% unique() %>% sort()
    cells <- depmap_ccle %>% dplyr::filter(tolower(primary_disease) %in% tolower(ccle_cancer_types))

    # cancer subtypes
    if (subtypes_optional_default_all) {
      ccle_cancer_subtypes <- subtypes
    } else {
      subtype_not_found <- list()
      subtype_found <- list()
      for (subtype in ccle_cancer_subtypes) {
        subtype <- str_trim(subtype)
        if (tolower(subtype) %in% tolower(subtypes)) {
          subtype_found <- append(subtype_found, subtype)
        } else {
          subtype_not_found <- append(subtype_not_found, subtype)
        }
      }
      if (length(subtype_found) == 0) {
        return (error_response(res, project_name = "DepMap", error_type="no_queried_cancer_subtypes_found"))
      } else {
        ccle_cancer_subtypes <- subtype_found
        # if there is only one subtype selected then apply filter on cells
        if(length(ccle_cancer_subtypes)==1){
            cells <- cells %>%dplyr::filter(tolower(Subtype) %in% tolower(ccle_cancer_subtypes))
        }
      }
    }
    # if there are valid subtypes selected (all or multiple or one)
    cells_names <- cells[["CCLE_Name"]]
    cells <- cells[["patient_id"]]
    names(cells) <- cells_names

    # CCLE cell lines
    if (cell_lines_optional_default_all) {
      ccle_cells <- cells
    } else {
      cell_line_not_found <- list()
      cell_line_found <- list()
      for (c in cell_line) {
        c <- str_trim(c)
        if (tolower(c) %in% tolower(cells)) {
          cell_line_found <- append(cell_line_found, c)
        } else {
          cell_line_not_found <- append(cell_line_not_found, c)
        }
      }
      if (length(cell_line_found) == 0) {
        return (error_response(res, project_name = "DepMap", error_type="no_queried_cell_lines_found"))
      } else {
        ccle_cells <- cell_line_found
      }
    }
    # print(length(ccle_cells))

    # DepMap gene/drug selection
    if (tolower(gene1) %in% depmap_genes) {
      print("test")
    } else {
      return (error_response(res, project_name = "DepMap", error_type="queried_gene_not_found"))
    }
  }

  list(analysis_parameter = paste0("The message is: '", cbind(unlist(values(parameter_dict))), "'"))
}


# #* only get data on list
# #* @filter data_check
# function(req){
#   cat(as.character(Sys.time()), "-",
#       req$QUERY_STRING, "\n")
#   plumber::forward()
# }


# #* Get the data inside of project_data
# #* @param file:[str] subproject name  + file name
# #* @serializer contentType list(type="text/csv")
# #* @get /project_data/<file>
# function(file = "DepMap_Expression", res){
#   #Set work directory to current file
#   setwd("/Applications/Codes/cSurvival_dev")
#   #!!! subject to change
#   #split the query into two parts: subproject(folder)
#   subproject = unlist(str_split(file, "_", n = 2))[1]
#   #and filename(to get the exact file based on query command)
#   filename = unlist(str_split(file, "_", n = 2))[2]
#   #Special case for mutation, difference between DepMap and Other
#   if(filename == "Mutation" & subproject == "DepMap"){
#     filename = command_dictionary[["DepMap_Mutation"]]
#   }else{
#     filename = command_dictionary[[filename]]  
#   }
#   filename <- paste0("./www/project_data/",subproject,"/",filename)
#   cat("The directory api is trying to access is: ",filename,"\n")
#   plumber::include_file(filename, res, "text/csv")
# }
# 
# # unlist(strsplit("DepMap+ccle_df.csv","+", fixed = TRUE))[2]
# 
# 


# #* Get the data inside of project_data
# #* @param subproject:[str] subproject folder name
# #* @param file:[str] file name
# #* @serializer contentType list(type="text/csv")
# #* @get /project_data/<subproject>/<file>
# function(subproject = "DepMap",file = "ccle.csv", res){
#   #Set work directory to current file
#   setwd("/Applications/Codes/cSurvival_dev")
#   #!!! subject to change
#   filename <- paste0("./www/project_data/",subproject,"/",file)
#   # data <- read.csv(filename)
#   res$setHeader("foo", "bar")
#   plumber::include_file(filename, res, "text/csv")
#   #plumber::as_attachment(data, paste0(subproject,"_",file))
# }
