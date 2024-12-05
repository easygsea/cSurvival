current_version = "v1.0.6"
d_red = "#EA5455"
bcols = c("#f4f4f4","#f5f0ed") # colors for parameter wellpanels EFF5FA EAEDF3
#the maximum -log10 value percentile heatmaps will apply as darkest color
heatmap_maximum_thershold = 3

# ------ available projects -----
pro_dir <- paste0(getwd(),"/www/project_data/")
projects <- read_csv(paste0(pro_dir,"project_name.csv")) %>% arrange(project_id)

projects_abbr <- projects$project_id
projects_full <- projects$name
names(projects_abbr) <- paste0(projects_abbr,": ", projects_full)

projects_cat <- str_split(projects_abbr,"-") %>% lapply(., function(x) x[1]) %>% unlist()
projects_cats <- unique(projects_cat)

projects <- projects_cats %>% lapply(., function(x){
  id <- x == projects_cat; prs <- projects_abbr[id]; return(prs) })
names(projects) <- projects_cats

projects[["DepMap"]][["DepMap-CRISPR: Chronos CRISPR Data"]] <- "DepMap-CRISPR"
projects[["DepMap"]][["DepMap-Drug: Drug sensitivity"]] <- "DepMap-Drug"
projects[["DepMap"]][["DepMap-RNAi: RNAi screening"]] <- "DepMap-RNAi"

projects <- sort_list(projects)

# ------- CCLE basic info --------
df_ccle <- fread(paste0(pro_dir,"DepMap/ccle.csv"), select = c("patient_id","CCLE_Name","gender","primary_or_metastasis","primary_disease","Subtype"))

