bcols = c("#FFFDE7","#FFEBEE") # colors for parameter wellpanels

# # available GMT gene sets
# gmt_dir <- paste0(dirname(getwd()),"/eVITTA/easyGSEA/www/gmts/hsa")
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