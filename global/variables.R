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
