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
