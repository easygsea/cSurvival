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



wget("https://ndownloader.figshare.com/files/26261527", '/Applications/Codes/cSurvival/project_data/DepMap/mutation.csv')
## create ExperimentHub query object
eh <- ExperimentHub()
query(eh, "depmap")

## retrieve metadata about cancer cell lines
metadata <- depmap::depmap_metadata()
print(metadata)




#FUNCTION----
#takes in a processed df and return the complete list of 
generate_full_name_table <- function(df){
   egENS <- toTable(org.Hs.egENSEMBL)
   egSYMBOL <- toTable(org.Hs.egSYMBOL)
   # bind the tables
   id_table <- egENS %>% left_join(egSYMBOL, by = "gene_id")
   
   df_columns <- colnames(df)
   for(index in 1:length(df_columns)){
      row_number <- which(id_table$symbol %in% df_columns[index])
   #If there is no matching, just continue
      if (!length(row_number)){
         next
      }
      temp_gene_id <- id_table$gene_id[row_number]
      temp_ensem <- id_table$ensembl_id[row_number]
   
      temp_column <- df_columns[index]
      temp_column <- paste0(temp_column,'|',unique(temp_gene_id))
      for(ele_ensem in temp_ensem){
         temp_column <- paste0(temp_column,'|',ele_ensem)
      }
      df_columns[index] <- temp_column
   }
   egENS <- NULL
   egSYMBOL <- NULL
   id_table <- NULL
   return(df_columns)
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
         dplyr::select(patient_id, everything())
   )
   if(!inherits(df_mir_scale, "try-error")){
      unlink(df_mir_scale_path, recursive = T)
   }
   # output the df to directory
   fwrite(df_mir_scale, file = df_mir_scale_path)
}


#START OF DF_GENE----
#TPM <- eh[["EH2264"]]
TPM <- depmap::depmap_TPM()
#We want to decrease the size of TPM
TPM <- TPM[, c("depmap_id", "gene",'expression')]


aftersplit <- split(TPM,unique(TPM$depmap_id))

#for(index in 1:length(aftersplit)){
for(index in 1:length(aftersplit)){
   temp <- data.frame(aftersplit[index])
   colnames(temp) <- c("depmap_id", "gene",'expression')
   print(index)
   name_list <- colnames(temp)
   new_col_name <- c('depmap_id',temp$gene)
   
   if(index == 1){
     df_gene_total <- reshape(temp, idvar = name_list[1], timevar = name_list[2], direction = "wide")
     colnames(df_gene_total) <- new_col_name
   }
   else{
      df_gene <- reshape(temp, idvar = name_list[1], timevar = name_list[2], direction = "wide")
      colnames(df_gene) <- new_col_name
      df_gene_total <- rbind(df_gene_total,df_gene)
   }
}
TPM <- NULL

#START OF Computer Specific because I lost my environment during my session
write.csv(df_gene_total,'df_gene.csv')
df_gene_total <- read.csv("/Applications/Codes/cSurvival/basic_scripts/df_gene.csv")
df_gene_total <- df_gene_total[ , -which(names(df_gene_total) %in% c("X"))]
change_name_list = colnames(df_gene_total)
#Split the column names
for (index in 1 : length(change_name_list)){
   change_name_list[index] <- strsplit(change_name_list[index],split = '.',fixed = TRUE)[[1]][1]
}

#Change column names of df_gene_total
colnames(df_gene_total) <- change_name_list
colnames(df_gene_total)[1] <- 'patient_id'
colnames(df_gene_total) <- generate_full_name_table(df_gene_total)


write.csv(df_gene_total,'df_gene.csv',row.names = FALSE)


#write scale df
df_gene_path <- "/Applications/Codes/cSurvival/project_data/DepMap/df_gene.csv"
df_gene_scale_path <- "/Applications/Codes/cSurvival/project_data/DepMap/df_gene_scale.csv"

scale_a_df(df_gene_path,df_gene_scale_path)

rm(df_gene_path)
rm(df_gene_scale_path)

#TEST WORKS
#test <- read.csv('df_gene.csv')
#END OF Computer Specific



#START OF DF_SNV----
#file path
filename = "/Applications/Codes/cSurvival/basic_scripts/mutation.csv"
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
      temp <- paste0(temp, '|', temp_mutation)
      #print(temp)
      df_snv_class[row_number,temp_gene] <- temp
   }
   if(is.na(df_snv_type[row_number,temp_gene])){
      df_snv_type[row_number,temp_gene] <- temp_type
   }
   else{
      #get current value of the cell
      temp <- df_snv_type[row_number,temp_gene]
      temp <- paste0(temp, '|', temp_type)
      #print(temp)
      df_snv_type[row_number,temp_gene] <- temp
   }
}

colnames(df_snv_class) <- generate_full_name_table(df_snv_class)
colnames(df_snv_type) <- generate_full_name_table(df_snv_type)

#df ready, now write csv
write.csv(df_snv_class,'df_snv_class.csv',row.names = FALSE)
write.csv(df_snv_type,'df_snv_type.csv',row.names = FALSE)

rm(df_snv_class)
rm(df_snv_type)
rm(df_mutations)


#START OF DF_CNV----
raw_data <- eh[["EH2262"]]#depmap::depmap_copyNumber()

raw_data <- raw_data[, c("depmap_id", "gene_name", "log_copy_number")]
raw_data <- data.frame(raw_data)

patients <- unique(raw_data$depmap_id)
change_name_list <- unique(raw_data$gene_name)
change_name_list <- str_sort(change_name_list, numeric = TRUE)

#The df we are going to write as csv
df_cnv <- data.frame(matrix(NA, nrow = length(patients), ncol = 1+length(change_name_list)))

colnames(df_cnv) <- c('patient_id',change_name_list)
df_cnv$patient_id <- patients
#Initialize the same empty type df
rm(patients)


#do first 5 million because last time I ran the server blacked out
for(index in 1: 5000000){
   print(index)
   temp_gene <- raw_data[index,"gene_name"]
   temp_patient <- raw_data[index,"depmap_id"]
   temp_num <- raw_data[index,"log_copy_number"]
   row_number <- which(df_cnv$patient_id %in% temp_patient)
   df_cnv[row_number,temp_gene] <- temp_num
}

#5 million to 10 million   
for(index in 5000001: 10000000){
   print(index)
   temp_gene <- raw_data[index,"gene_name"]
   temp_patient <- raw_data[index,"depmap_id"]
   temp_num <- raw_data[index,"log_copy_number"]
   row_number <- which(df_cnv$patient_id %in% temp_patient)
   df_cnv[row_number,temp_gene] <- temp_num
}



rm(raw_data)
rm(temp_gene)
rm(temp_num)
rm(temp_patient)



# aftersplit <- split(TPM,unique(TPM$depmap_id))

#for(index in 1:length(aftersplit)){
# for(index in 1:length(aftersplit)){
#    temp <- data.frame(aftersplit[index])
#    colnames(temp) <- c("depmap_id", "gene",'expression')
#    print(index)
#    name_list <- colnames(temp)
#    new_col_name <- c('depmap_id',temp$gene)
#    
#    if(index == 1){
#       df_gene_total <- reshape(temp, idvar = name_list[1], timevar = name_list[2], direction = "wide")
#       colnames(df_gene_total) <- new_col_name
#    }
#    else{
#       df_gene <- reshape(temp, idvar = name_list[1], timevar = name_list[2], direction = "wide")
#       colnames(df_gene) <- new_col_name
#       df_gene_total <- rbind(df_gene_total,df_gene)
#    }
# }
# TPM <- NULL
