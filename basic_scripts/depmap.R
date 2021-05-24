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
#FUNCTION----
#takes in a processed df and return the complete list of 
generate_full_name_table <- function(df, ENSG = TRUE){
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
      if(ENSG == TRUE){
         for(ele_ensem in temp_ensem){
            temp_column <- paste0(temp_column,'|',ele_ensem)
         }
      }
      
      df_columns[index] <- temp_column
   }
   egENS <- NULL
   egSYMBOL <- NULL
   id_table <- NULL
   return(df_columns)
}



generate_full_name_table_SNV <- function(df){
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
      temp_column <- df_columns[index]
      temp_column <- paste0(temp_column,'|',unique(temp_gene_id))
      df_columns[index] <- temp_column
   }
   
   
   egENS <- NULL
   egSYMBOL <- NULL
   id_table <- NULL
   return(df_columns)
}


#generate full name for gene_name column
generate_full_name_column_crispr <- function(df){
   egENS <- toTable(org.Hs.egENSEMBL)
   egSYMBOL <- toTable(org.Hs.egSYMBOL)
   # bind the tables
   id_table <- egENS %>% left_join(egSYMBOL, by = "gene_id")
   
   df_columns <- unique(df$gene_name)
   for(index in 1:length(df_columns)){
      print(index)
      row_number <- which(id_table$symbol %in% df_columns[index])
      #If there is no matching, just continue
      if (!length(row_number)){
         next
      }
      temp_gene_id <- id_table$gene_id[row_number]
      #temp_ensem <- id_table$ensembl_id[row_number]
      
      temp_column <- df_columns[index]
      temp_column <- paste0(temp_column,'|',unique(temp_gene_id))
      # for(ele_ensem in temp_ensem){
      #   temp_column <- paste0(temp_column,'|',ele_ensem)
      # }
      df_columns[index] <- temp_column
      print(temp_column)
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


#This function takes in a df and returns a list
#that stores different indices where the genes differ
find_start_of_genes <- function(df,mode = 'crispr'){
   #Use a variable called target to measure the index
   if(mode == 'drug_sensitivity'){
      target = 'compound'
   }
   if(mode == 'crispr'){
      target = 'gene_name'
   }
   for(index in 1:nrow(df)){
      current <- df[index,target]
      #Initialize gene so far and result list   
      if(index == 1){
         gene_so_far = df[index,target]
         result <- c(index)
      }
      #current still matches gene so far
      if(current != gene_so_far){
         #update gene so far
         gene_so_far = current
         result <- c(result,index)
      }
      
   }
   return(result)
}
#Returns a list of gender for df, which is based on gender_df
add_gender <- function(df,gender_df){
   gender_list <- c()
   for(index in 1:nrow(df)){
      row_number <- which(gender_df$depmap_id %in% df[index,1])
      if(!length(row_number)){
         gender_list[index] <- 'Unknown'
      }
      else{
         gender_list[index] <- gender_df$sex[row_number]
      }
      print(index)
   }
   
   return(gender_list)
}



################### Preteomic Specific Function

#This function is used to split the patient_id part of df-proteomic 
#and we can use the result to look up in CCLE.csv
split_TenPx <- function(df){
   result <- strsplit(df$patient_id,split = '_TenPx')
   for(index in 1: length(result)){
      result[index] = result[[index]][1]
   }
   result <- unlist(result)
   df$patient_id <- result
   return(df)
}

#Returns a updated proteomic_df where the patient_id completed by looking up ccle.csv
ccle_lookup <- function(df,df_ccle){
   df <- split_TenPx(df)
   result <- c()
   for(index in 1:nrow(df)){
      row_number <- which(df_ccle$CCLE_Name %in% df[index,1])
      if(!length(row_number)){
         #Try get rid of the X at the beginnning of the strig and try again
         df[index,1] <- substr(df[index,1], 2, nchar(df[index,1]))
         row_number <- which(df_ccle$CCLE_Name %in% df[index,1])
         #Still cannot find
         if(!length(row_number)){
            print(paste('Patient ID not found and index number is: ', index))
         }
         #Now can find
         else{
            result[index] <- df[index,1]
         }
         #print(paste('Patient ID not found and their ID is: ', df[index,1]))
         #if not found, keep original
         result[index] <- df[index,1]
      }
      else{
         result[index] <- df_ccle$patient_id[row_number]
      }
      
   }
   df$patient_id <- result
   return(df)
}




## create ExperimentHub query object
eh <- ExperimentHub()
query(eh, "depmap")

## retrieve metadata about cancer cell lines
metadata <- depmap::depmap_metadata()

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


#STATR OF DF_GENE----
TPM <- depmap::depmap_TPM()
#We want to decrease the size of TPM
raw_data <- TPM[, c("depmap_id", "gene",'expression')]
rm(TPM)
raw_data <- data.frame(raw_data)

patients <- unique(raw_data$depmap_id)
change_name_list <- unique(raw_data$gene_name)
change_name_list <- str_sort(change_name_list, numeric = TRUE)

#The df we are going to write as csv
dmp_df_gene <- data.frame(matrix(NA, nrow = length(patients), ncol = 1+length(change_name_list)))

colnames(dmp_df_gene) <- c('patient_id',change_name_list)
dmp_df_gene$patient_id <- patients
#Initialize the same empty type df
rm(patients)


#do first half first because last time I ran the server blacked out
for(index in 1: (1/2*nrow(raw_data))){
   print(index)
   temp_gene <- raw_data[index,"gene"]
   temp_patient <- raw_data[index,"depmap_id"]
   temp_num <- raw_data[index,"expression"]
   row_number <- which(df_gene$patient_id %in% temp_patient)
   dmp_df_gene[row_number,temp_gene] <- temp_num
}
#You can save this just in case to save our progress
write.csv(dmp_df_gene,'~/ShinyApps/www/project_data/DepMap/df_gene.csv',row.names = FALSE)
#pause and run
for(index in ((1/2*nrow(raw_data))+1): nrow(raw_data)){
   print(index)
   temp_gene <- raw_data[index,"gene"]
   temp_patient <- raw_data[index,"depmap_id"]
   temp_num <- raw_data[index,"expression"]
   row_number <- which(df_gene$patient_id %in% temp_patient)
   dmp_df_gene[row_number,temp_gene] <- temp_num
}
#Save as csv
write.csv(dmp_df_gene,'~/ShinyApps/www/project_data/DepMap/df_gene.csv',row.names = FALSE)
#remove variables in the environment
rm(temp_gene)
rm(temp_num)
rm(temp_patient)
rm(dmp_df_gene)


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


#do first 10 million because last time I ran the server blacked out

for(index in 1: 10000000){
   print(index)
   temp_gene <- raw_data[index,"gene_name"]
   temp_patient <- raw_data[index,"depmap_id"]
   temp_num <- raw_data[index,"log_copy_number"]
   row_number <- which(df_cnv$patient_id %in% temp_patient)
   df_cnv[row_number,temp_gene] <- temp_num
}


for(index in 10000000: 20000000){
   print(index)
   temp_gene <- raw_data[index,"gene_name"]
   temp_patient <- raw_data[index,"depmap_id"]
   temp_num <- raw_data[index,"log_copy_number"]
   row_number <- which(df_cnv$patient_id %in% temp_patient)
   df_cnv[row_number,temp_gene] <- temp_num
}


for(index in 20000000: 30000000){
   print(index)
   temp_gene <- raw_data[index,"gene_name"]
   temp_patient <- raw_data[index,"depmap_id"]
   temp_num <- raw_data[index,"log_copy_number"]
   row_number <- which(df_cnv$patient_id %in% temp_patient)
   df_cnv[row_number,temp_gene] <- temp_num
}

for(index in 30000000: nrow(raw_data)){
   print(index)
   temp_gene <- raw_data[index,"gene_name"]
   temp_patient <- raw_data[index,"depmap_id"]
   temp_num <- raw_data[index,"log_copy_number"]
   row_number <- which(df_cnv$patient_id %in% temp_patient)
   df_cnv[row_number,temp_gene] <- temp_num
}
View(df_cnv)

colnames(df_cnv) <- generate_full_name_table(df_cnv)
write.csv(df_cnv,'~/ShinyApps/www/project_data/DepMap/df_cnv.csv',row.names = FALSE)
View(df_cnv)
test <- read.csv("DepMap/df_cnv.csv")
View(test)


#START OF PROTEOMIC----
setwd("~/ShinyApps/basic_scripts")
#gets stored in basic_scripts folder after set working directory to this one
wget("https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz")
filename = "~/ShinyApps/basic_scripts/protein_quant_current_normalized.csv.gz"
# load to a dataframe
raw_data <- fread(file = filename,header=T)#,select = c('Hugo_Symbol', 'DepMap_ID','Variant_Classification','Variant_Type'))
raw_data <- data.frame(raw_data)
raw_data$gene_name <- paste0(raw_data$Gene_Symbol, "|", raw_data$Uniprot,'|', raw_data$Uniprot_Acc)

df_proteomic <- dplyr::select(raw_data, matches("_TenPx"))
df_proteomic <- cbind(gene_name = raw_data$gene_name, df_proteomic)

rownames(df_proteomic) <- df_proteomic$gene_name
df_proteomic <- as.data.frame(t(as.matrix(df_proteomic)))
df_proteomic$patient_id = rownames(df_proteomic)

rownames(df_proteomic) <- NULL
df_proteomic <- df_proteomic %>% dplyr::select('patient_id', everything())
df_proteomic <- df_proteomic[-1,]

#Get the ccle look up so that we can have ACH- patient id for patient id column
df_ccle <- fread(file = '~/ShinyApps/www/project_data/DepMap/ccle.csv', header = T, select = c('CCLE_Name','patient_id'))
df_proteomic <- ccle_lookup(df = df_proteomic, df_ccle = df_ccle)

#Now we found that there are some space in the values and we want to get rid of them
#substituting the values using gsub()
for(col in colnames(df_proteomic)){
   df_proteomic[,col] <- gsub(' ','',df_proteomic[,col], fixed = TRUE)
}


#df_proteomic$`SLC12A2|S12A2_HUMAN|P55011` <- gsub(' ','',df_proteomic$`SLC12A2|S12A2_HUMAN|P55011`, fixed = TRUE)
fwrite(df_proteomic, file = '~/ShinyApps/www/project_data/DepMap/df_proteomic.csv')
rm(df_proteomic)
rm(raw_data)
rm(df_ccle)







#START OF CRISPR-Cas9----
raw_data <- depmap::depmap_crispr()
raw_data <- raw_data[, c("depmap_id", "gene_name", "dependency")]
raw_data <- data.frame(raw_data)
#dropped na
raw_data <- raw_data[!(is.na(raw_data$dependency)), ]
#Order the df by gene name, and we hope to solve this question just in one time
raw_data <- raw_data[order(raw_data$gene_name),]
#Add gender
gender <- depmap::depmap_metadata()
gender <- gender[, c('depmap_id','sex')]


gender_list <- add_gender(df = raw_data,gender_df = gender)
gender_list

colnames(raw_data)[1] <- 'patient_id'

#This takes a little while
#This list will be used to know which rows to slice
slice_list <- find_start_of_genes(raw_data)
raw_data <- cbind(raw_data, gender=gender_list)
colnames(raw_data)[4] = 'gender'

#Double check and the slice list matches
#length(slice_list)
#length(unique(raw_data$gene_name))
#First get full gene names sorted out:

full_gene_name <- generate_full_name_column_crispr(raw_data)
print(length(full_gene_name))
print(length(slice_list))
full_gene_name

#Create directory
for(index in 1:length(full_gene_name)){
   temp_path = paste0('~/ShinyApps/www/project_data/DepMap/CRISPR-Cas9/',full_gene_name[index],'/')
   dir.create(path = temp_path)
}


#This part takes in a df and a index list and writes in different folders
#have to separate into many parts because of massive runtime


for(index in 1:length(slice_list)){
   #(1/3*length(slice_list))){
   # if(index == 1){
   #   start_index <- 1
   #   end_index <- (slice_list[2] - 1)
   # }
   #If reaches the end
   print(index)
   if(index == length(slice_list)){
      start_index <- slice_list[index]
      end_index <- nrow(raw_data)
   }
   else{
      start_index <- slice_list[index]
      end_index <- slice_list[index + 1]
   }
   temp_df <- raw_data[start_index:end_index,]
   temp_df$gene_name <- full_gene_name[index]
   temp_path = paste0('~/ShinyApps/www/project_data/DepMap/CRISPR-Cas9/',full_gene_name[index],'/','df_survival.csv')
   fwrite(temp_df, file = temp_path)
   #write.csv(temp_df, temp_path, row.names = FALSE)
}

map_path = '~/ShinyApps/www/project_data/DepMap/CRISPR-Cas9/directory.csv'
full_gene_name <- data.frame(full_gene_name)
colnames(full_gene_name) <- 'gene_name'
fwrite(full_gene_name, file = map_path)
rm(map_path)
rm(raw_data)
rm(slice_list)
rm(full_gene_name)
rm(gender)
rm(gender_list)


#START OF Drug_sensitivity----
raw_data <- depmap::depmap_drug_sensitivity()
raw_data <- raw_data[, c("depmap_id", "compound",'dependency')]
raw_data <- data.frame(raw_data)
raw_data <- raw_data[!(is.na(raw_data$dependency)), ]

raw_data <- raw_data[order(raw_data$compound),]
#Add gender
gender <- depmap::depmap_metadata()
gender <- gender[, c('depmap_id','sex')]


gender_list <- add_gender(df = raw_data,gender_df = gender)
gender_list

colnames(raw_data)[1] <- 'patient_id'

#This takes a little while
#This list will be used to know which rows to slice
slice_list <- find_start_of_genes(raw_data,mode = 'drug_sensitivity')
raw_data <- cbind(raw_data, gender=gender_list)
rm(gender)
rm(gender_list)

#Double check and the slice list matches
#length(slice_list)
#length(unique(raw_data$gene_name))
#First get full gene names sorted out:

full_gene_name <- unique(raw_data$compound)
print(length(full_gene_name))
print(length(slice_list))
full_gene_name

#Create directory
for(index in 1:length(full_gene_name)){
   temp_path = paste0('~/ShinyApps/www/project_data/DepMap/Drug_sensitivity/',full_gene_name[index],'/')
   dir.create(path = temp_path)
}


#This part takes in a df and a index list and writes in different folders
#have to separate into many parts because of massive runtime


for(index in 1:length(slice_list)){
   #If reaches the end
   print(index)
   if(index == length(slice_list)){
      start_index <- slice_list[index]
      end_index <- nrow(raw_data)
   }
   else{
      start_index <- slice_list[index]
      end_index <- slice_list[index + 1]
   }
   temp_df <- raw_data[start_index:end_index,]
   temp_path = paste0('~/ShinyApps/www/project_data/DepMap/Drug_sensitivity/',full_gene_name[index],'/','df_survival.csv')
   fwrite(temp_df, file = temp_path)
   #write.csv(temp_df, temp_path, row.names = FALSE)
}

map_path = '~/ShinyApps/www/project_data/DepMap/Drug_sensitivity/directory.csv'
full_gene_name <- data.frame(full_gene_name)
colnames(full_gene_name) <- 'compound'
fwrite(full_gene_name, file = map_path)
rm(map_path)
rm(raw_data)
rm(slice_list)
rm(full_gene_name)


#START OF RNAi----
raw_data <- depmap::depmap_rnai()
raw_data <- raw_data[, c("depmap_id", 'gene_name','dependency')]
raw_data <- data.frame(raw_data)
raw_data <- raw_data[!(is.na(raw_data$dependency)), ]
raw_data <- raw_data[order(raw_data$gene_name),]


gender <- depmap::depmap_metadata()
gender <- gender[, c('depmap_id','sex')]



#gender_list <- c()
#gender <- data.frame(gender)




gender_list <- add_gender(df = raw_data,gender_df = gender)
gender_list

colnames(raw_data)[1] <- 'patient_id'

#This takes a little while
#This list will be used to know which rows to slice
slice_list <- find_start_of_genes(raw_data)
raw_data <- cbind(raw_data, gender=gender_list)
colnames(raw_data)[4] = 'gender'

#Double check and the slice list matches
#length(slice_list)
#length(unique(raw_data$gene_name))
#First get full gene names sorted out:

full_gene_name <- generate_full_name_column_crispr(raw_data)
print(length(full_gene_name))
print(length(slice_list))
full_gene_name

#Create directory
for(index in 1:length(full_gene_name)){
   temp_path = paste0('~/ShinyApps/www/project_data/DepMap/RNAi/',full_gene_name[index],'/')
   dir.create(path = temp_path)
}


#This part takes in a df and a index list and writes in different folders
#have to separate into many parts because of massive runtime


for(index in 1:length(slice_list)){
   #(1/3*length(slice_list))){
   # if(index == 1){
   #   start_index <- 1
   #   end_index <- (slice_list[2] - 1)
   # }
   #If reaches the end
   print(index)
   if(index == length(slice_list)){
      start_index <- slice_list[index]
      end_index <- nrow(raw_data)
   }
   else{
      start_index <- slice_list[index]
      end_index <- slice_list[index + 1]
   }
   temp_df <- raw_data[start_index:end_index,]
   temp_df$gene_name <- full_gene_name[index]
   temp_path = paste0('~/ShinyApps/www/project_data/DepMap/RNAi/',full_gene_name[index],'/','df_survival.csv')
   fwrite(temp_df, file = temp_path)
   #write.csv(temp_df, temp_path, row.names = FALSE)
}

map_path = '~/ShinyApps/www/project_data/DepMap/RNAi/directory.csv'
full_gene_name <- data.frame(full_gene_name)
colnames(full_gene_name) <- 'gene_name'
fwrite(full_gene_name, file = map_path)
rm(map_path)
rm(raw_data)
rm(slice_list)
rm(full_gene_name)
rm(gender)
rm(gender_list)

