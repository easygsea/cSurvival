# directory to temporarily save variables
surv_dir <- paste0(dirname(dirname(getwd())), "/variables/")

# ------ available GMT gene sets -----
# gmt_dir <- paste0(dirname(getwd()),"/eVITTA_dev/easyGSEA/www/gmts/hsa/")
# gmt_dir <- "http://tau.cmmt.ubc.ca/eVITTA/easyGSEA/gmts/hsa/"
gmt_dir <- "/home/eVITTA/ShinyApps/easyGSEA/www/gmts/hsa/" # tau path
# # collection in a df
# in_gmt_lst <- paste0(dirname(getwd()),"/eVITTA_dev/easyGSEA/www/gmts/gmts_list.csv") # code test locally
# in_gmt_lst <- "http://tau.cmmt.ubc.ca/eVITTA/easyGSEA/gmts/gmts_list.csv"
in_gmt_lst <- "/home/eVITTA/ShinyApps/easyGSEA/www/gmts/gmts_list.csv" # tau path
gmt_collection_df <- read_csv(in_gmt_lst,col_names = F) %>% dplyr::filter(X1 == "hsa")
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
  paste0(gmt_dir,sub_dir,"/",db) %>% gsub(" ","%20",.)
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
# ------- names for input modes ----------
input_mode_names <- c(
  "Expression (FPKM)" = "rna"
  ,"Mutation" = "snv"
  ,"Copy number" = "cnv"
  ,"Expression (RPM)" = "mir"
  ,"Methylation beta-value" = "met"
  ,"Mean Z-score" = "lib"
  ,"Mean Z-score" = "manual"
)

# ------- survival analysis methods ---------
surv_methods <- c("Cox proportional-hazards (PH) model"="cox", "Kaplan-Meier (KM) log rank"="km")

# ------- gene set data types, library or manual -------
data_types_gs <- c("Gene set"="lib",
                   "Gene set (manual)"="manual")

# ------- placeholder for gene search field ----------
g_placeholder <- "On top left, select a project(s) and click the confirmation button to load genes ..."

# ------- algorithms for somatic variant calling, e.g. mutect -----------
snv_algorithms <- list(
  "MuTect" = "mutect"
  ,"VarScan" = "varscan"
  ,"SomaticSniper" = "somaticsniper"
  ,"MuSE" = "muse"
)

# ------ variant types/classes, for classification into non- or synonymous mutations -------
variant_types <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins"
                  ,"Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"
                  ,"Splice_Site","Translation_Start_Site"
                  ,"Silent","Splice_Region","Intron","5'UTR","RNA","3'UTR"        
                  ,"5'Flank","3'Flank","IGR")

variant_types_non <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins"
                       ,"Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"
                       ,"Splice_Site","Translation_Start_Site")

variant_types_syn <- c("Silent","Splice_Region","Intron","5'UTR","RNA","3'UTR"        
                       ,"5'Flank","3'Flank","IGR")

# ------- dynamic variables need initiatialization and updates according to inputs' changes ------
dyn_list <- function(x){
  list(
    # category to analyze
    cat_id <- paste0("cat_",x)
    # type of db to analyze
    ,db_id <- paste0("db_",x)
    # gene to analyze
    ,g_ui_id <- paste0("g_",x)
    # gene set to analyze
    ,gs_mode_id <- paste0("gs_mode_",x)
    ,gs_db_id <- paste0("gs_db_",x)
    ,gs_lib_id <- paste0("gs_l_",x)
    ,gs_lib_genes_id <- paste0("gs_lgs_",x) # verbtxt
    # user input gene(s) to filter gs
    ,gs_gene_id <- paste0("gs_lg_",x)
    # # search btn click
    # ,gs_gene_btn_id <- paste0("gs_lg_",x,"_search")
    ,gs_gene_genes_id <- paste0("gs_lgg_",x) # verbtxt
    # manual gene input
    ,gs_manual_id <- paste0("gs_m_",x)
    ,gs_manual_btn <- paste0("add_btn_",x)
    ,gs_genes_id <- paste0("gs_mg_",x)
    ,lower_id <- paste0("lower_",x)
    ,higher_id <- paste0("upper_",x)
    ,step_id <- paste0("step_",x)
    ,snv_id <- paste0("snv_method_",x)
    ,non_id <- paste0("nonsynonymous_",x)
    # ,syn_id <- paste0("synonymous_",x)
    ,iter_id <- paste0("iter_",x)
    ,clow_id <- paste0("clow_",x)
    ,cnv_id <- paste0("cnv_par_",x)
  )
}
# the data that tell what Target projects data have
TARGET_existing_data <- fread(paste0(getwd(),"/project_data/TARGET_existing_data.csv"), sep = ",") %>%
  mutate(existing_data = str_split(existing_data, pattern = ","))

# -------- TCGA survival endpoints ----------
tcga_stypes <- c(
  "Overall survival (OS)" = "OS"
  ,"Progression-free interval (PFI)" = "PFI"
  ,"Disease-free interval (DFI)" = "DFI"
  # ,"Progression-specific survival (PSS)" = "pss"
  ,"Disease-specific survival (DSS)" = "DSS"
)

# -------- TCGA clinical codes ----------
tcga_codes <- fread(paste0(getwd(),"/project_data/tcga_codes.tsv"), sep = "\t")

codes_color <- list(
  "yes" = "#339900" #99cc33
  ,"no" = "#cc3300"
  ,"caution" = "#f58c00" # ff9966
  ,"app|caution" = "#f58c00"
  ,"app" = "#339900"
  ,"acc" = "#339900"
)

codes_icon <- list(
  "yes" = icon("check")
  ,"no" = icon("times")
  ,"caution" = icon("exclamation")
  ,"app|caution" = icon("exclamation")
  ,"app" = icon("check")
  ,"acc" = icon("check")
)

# ----- YMD time units -----
ymd_names <- c(
  "Days"="d",
  "Months"="m",
  "Years"="y"
)
ymd_unit <- c(
  "y" = 365.25
  ,"m" = 30.4375
  ,"d" = 1
)