library("depmap")
library("ExperimentHub")
library('HelpersMG')
#This is needed for gene name conversion
library("org.Hs.eg.db")
library(data.table)
library(stringr)
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
library(TCGAutils)
#START OF SNV----
#IMPORTANT:set working directory to this file
#Get csv needed for snv
setwd("~/ShinyApps/basic_scripts")
wget("https://ndownloader.figshare.com/files/26261527")

#I set working directory to this file
filename = "~/ShinyApps/basic_scripts/26261527"
# load to a dataframe
df_mutations <- fread(file = filename,header=T,select = c('Hugo_Symbol', 'DepMap_ID','Variant_Classification','Variant_Type'))
df_mutations <- data.frame(df_mutations)

patients <- unique(df_mutations$DepMap_ID)
change_name_list <- unique(df_mutations$Hugo_Symbol)
change_name_list <- str_sort(change_name_list, numeric = TRUE)

#The df we are going to write as csv
df_snv_class <- data.frame(matrix(NA, nrow = length(patients), ncol = 1+length(change_name_list)))

#include a more column for patient ID
colnames(df_snv_class) <- c('patient_id',change_name_list)
df_snv_class$patient_id <- patients
#Initialize the same empty type df
df_snv_type <- df_snv_class
rm(patients)



for(index in 1: nrow(df_mutations)){
  print(index)
  temp_gene <- df_mutations[index,'Hugo_Symbol']
  temp_patient <- df_mutations[index,'DepMap_ID']
  
  temp_mutation <- df_mutations[index,'Variant_Classification']
  temp_type <- df_mutations[index,'Variant_Type']
  #Because the type df looks the same, so it shares the same row and column index
  row_number <- which(df_snv_class$patient_id %in% temp_patient)
  
  if(is.na(df_snv_class[row_number,temp_gene])){
    df_snv_class[row_number,temp_gene] <- temp_mutation
  }
  else{
    #get current value of the cell
    temp <- df_snv_class[row_number,temp_gene]
    #Add a check if already contains, then no need to append
    if(!(grepl(temp, temp_mutation))){
      temp <- paste0(temp, '|', temp_mutation)
    }
    #print(temp)
    df_snv_class[row_number,temp_gene] <- temp
  }
  if(is.na(df_snv_type[row_number,temp_gene])){
    df_snv_type[row_number,temp_gene] <- temp_type
  }
  else{
    #get current value of the cell
    temp <- df_snv_type[row_number,temp_gene]
    if(!(grepl(temp, temp_type))){
      temp <- paste0(temp, '|', temp_type)
    }
    #print(temp)
    df_snv_type[row_number,temp_gene] <- temp
  }
}

colnames(df_snv_class) <- generate_full_name_table(df_snv_class,ENSG = FALSE)
colnames(df_snv_type) <- generate_full_name_table(df_snv_type, ENSG = FALSE)

#df ready, now write csv
fwrite(df_snv_class, file = '~/ShinyApps/www/project_data/DepMap/df_snv_class.csv')
fwrite(df_snv_type, file = '~/ShinyApps/www/project_data/DepMap/df_snv_type.csv')
# write.csv(df_snv_class,'~/ShinyApps/www/project_data/DepMap/df_snv_class.csv',row.names = FALSE)
# write.csv(df_snv_type,'~/ShinyApps/www/project_data/DepMap/df_snv_type.csv',row.names = FALSE)

rm(df_snv_class)
rm(df_snv_type)
rm(df_mutations)