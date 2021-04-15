bcols = c("#f0eee9","#e9ebf0") # colors for parameter wellpanels

# ------ available projects -----
projects <- read_csv(paste0(getwd(),"/project_data/project_name.csv")) %>% arrange(project_id)

projects_abbr <- projects$project_id
projects_full <- projects$name
names(projects_abbr) <- paste0(projects_abbr,": ", projects_full)

projects_cat <- str_split(projects_abbr,"-") %>% lapply(., function(x) x[1]) %>% unlist()
projects_cats <- unique(projects_cat)

projects <- projects_cats %>% lapply(., function(x){
  id <- x == projects_cat; prs <- projects_abbr[id]; return(prs) })
names(projects) <- projects_cats

# ------ available GMT gene sets -----
gmt_dir <- paste0(dirname(getwd()),"/eVITTA_dev/easyGSEA/www/gmts/hsa/")
# collection in a df
gmt_collection_df <- read_csv(paste0(dirname(getwd()),"/eVITTA_dev/easyGSEA/www/gmts/gmts_list.csv"),col_names = F) %>% dplyr::filter(X1 == "hsa")
# available databases
gmt_dbs <- str_split(gmt_collection_df$X3,";") %>% lapply(function(x) gsub("_"," ",gsub("?[.]gmt$","",x)))
db_names <- gsub("_"," ",gsub("^\\d+_?","",gmt_collection_df$X2))
names(gmt_dbs) <- db_names
# function to find the full path to a selected GMT
retrieve_gmt_path <- function(db){
  db <- paste0(gsub(" ","_",db),".gmt")
  # determine which subdirectory the library is in
  i <- grep(db, gmt_collection_df$X3)
  sub_dir <- gmt_collection_df$X2[i]
  paste0(gmt_dir,sub_dir,"/",db)
}
# gmt_files <- list.files(gmt_dir, pattern = "[.]gmt$", recursive = T, full.names = T) %>% .[!grepl("TF2DNA[.]gmt$",.)]
# gmts <- gmt_files %>%
#   lapply(., function(x){
#     # cat_name <- basename(dirname(x)) %>% gsub("^\\d+_?","",.)
#     # db_name <- gsub("?[.]gmt$","",basename(x)) %>% as.list()
#     # names(db_name) <- cat_name
#     # return(db_name)
#     gmtPathways(x)
#   })
# names(gmts) <- gsub("?[.]gmt$","",basename(gmt_files))