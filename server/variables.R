# directory to temporarily save variables
surv_dir <- paste0(dirname(dirname(getwd())), "/variables/")

# ------ available GMT gene sets -----
# gmt_dir <- paste0(dirname(getwd()),"/eVITTA_dev/easyGSEA/www/gmts/hsa/")
gmt_dir <- "http://tau.cmmt.ubc.ca/eVITTA/easyGSEA/gmts/hsa/"
# gmt_dir <- "/home/eVITTA/ShinyApps/easyGSEA/www/gmts/hsa/" # tau path
# # collection in a df
# in_gmt_lst <- paste0(dirname(getwd()),"/eVITTA_dev/easyGSEA/www/gmts/gmts_list.csv") # code test locally
in_gmt_lst <- "http://tau.cmmt.ubc.ca/eVITTA/easyGSEA/gmts/gmts_list.csv"
# in_gmt_lst <- "/home/eVITTA/ShinyApps/easyGSEA/www/gmts/gmts_list.csv" # tau path
gmt_collection_df <- read_csv(in_gmt_lst,col_names = F) %>% dplyr::filter(X1 == "hsa")
# available databases
gmt_dbs <- str_split(gmt_collection_df$X3,";") %>% lapply(function(x) gsub("_"," ",gsub("?[.]gmt$","",x)))
db_names <- gsub("_"," ",gsub("^\\d+_?","",gmt_collection_df$X2))
names(gmt_dbs) <- db_names
gmt_dbs[["Other"]] <- gmt_dbs[["Other"]][-4]
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
# ------ flagged cases -----
flagged_cases <- fread(paste0(pro_dir,"/977/flagged_cases.tsv"),sep="\t") %>% .[["patient_id"]]
# ------- names for input modes ----------
input_mode_names <- c(
  "Expression (FPKM)" = "rna"
  ,"Mutation" = "snv"
  ,"Copy number" = "cnv"
  ,"Expression (RPM)" = "mir"
  ,"Methylation beta-value" = "met"
  ,"Mean Z-score" = "lib"
  ,"Mean Z-score" = "manual"
  ,"Median-centered RRPA value" = "rrpa"
  ,"Normalized protein expression" = "pro"
  ,"Gene effect (CERES)" = "crispr"
  ,"Gene effect (DEMETER2)" = "rnai"
  ,"Cell viability" = "drug"
)

# ------- survival analysis methods ---------
surv_methods <- c("Kaplan-Meier (KM) log rank"="km","Cox proportional-hazards (PH) model"="cox")

# ------- gene set data types, library or manual -------
data_types_gs <- c("Gene set"="lib",
                   "Gene set (manual)"="manual")

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
                  ,"5'Flank","3'Flank","IGR","WT")

variant_types_non <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins"
                       ,"Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"
                       ,"Splice_Site","Translation_Start_Site")

variant_types_syn <- c("Silent","Splice_Region","Intron","5'UTR","RNA","3'UTR"        
                       ,"5'Flank","3'Flank","IGR","WT")

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
    ,syn_id <- paste0("synonymous_",x)
    ,iter_id <- paste0("iter_",x)
    ,clow_id <- paste0("clow_",x)
    ,cnv_id <- paste0("cnv_par_",x)
    ,snv_uni_id <- paste0("snv_uni_",x)
  )
}
# ----- the data that tell what TCGA/TARGET has -----
TARGET_existing_data <- fread(paste0(pro_dir,"TARGET_existing_data.csv"), sep = ",") %>%
  mutate(existing_data = str_split(existing_data, pattern = ","))

TCGA_missing_data <- fread(paste0(pro_dir,"TCGA_missing_data.csv"), sep = ",") %>%
  mutate(missing_data = str_split(missing_data, pattern = ","))

# -------- TCGA survival endpoints ----------
tcga_stypes <- c(
  "Overall survival (OS)" = "OS"
  ,"Progression-free interval (PFI)" = "PFI"
  ,"Disease-free interval (DFI)" = "DFI"
  # ,"Progression-specific survival (PSS)" = "pss"
  ,"Disease-specific survival (DSS)" = "DSS"
)

# -------- TCGA clinical codes ----------
tcga_codes <- fread(paste0(pro_dir,"tcga_codes.tsv"), sep = "\t")

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

# -------- bstooltip help text --------
cat_id_q_txt <- "<b>Gene or locus</b>: To study if the expression level, mutational status, copy number, methylation, or protein level of a gene or locus correlates with poorer/better survival.<br><b>Gene set</b>: To study if the average expression level of a gene set correlates with cancer survival, e.g. genes in the same pathway, TF targets, drug targets, miRNA targets, interacting proteins, or user-defined list of genes."
db_id_q_txt <- "To study if cancer survival is associated with a gene\\'s expression level, mutational status, copy number variation; a microRNA\\'s expression; or the methylation level of a DNA segment."
g_id_txt <- paste0("Search and select. If a gene or locus is not found, try its alias names."
                   ," If still not found, it means its expression/alteration is barely detected in the selected cancer project."
                   ," Or, if proteomics, it has not been quantified."
                   ," Or, if pan-cancer analysis, its expression/alteration is not detected in all selected projects.")
cox_km_txt <- paste0("Select the method for analyzing and summarizing survival data. "
                     ,"KM describe the survival according to one factor under investigation; "
                     ,"Cox regression model assesses the effect of several risk factors simultaneously; "
                     ,"additional density and box plots are provided for DepMap data analysis and visualization."
)
padj_q_txt <- paste0("Adjustment for multiple testing arising from assessing a sequence of candidate thresholds with the minimum <i>P</i>-value method (Lausen and Schumacher, 1992)")

# ----- Miscellaneous ----
# correction methods for pairwise comparisons
pairwise_methods <- list(
  "Multiple comparisons test by Holm (1979)" = "holm"
  ,"Multiple comparisons test by Hochberg (1988)" = "hochberg"
  ,"Multiple comparisons test by Hommel (1988)" = "hommel"
  ,"Multiple comparisons test by Bonferroni correction" = "bonferroni"
  ,"Multiple comparisons test by Benjamini & Hochberg (1995)" = "BH"
  ,"Multiple comparisons test by Benjamini & Yekutieli (2001)" = "BY"
  ,"Multiple comparisons test by false discovery rate (FDR)" = "fdr"
)

# explanations for flagged cases
flagged_exp <- paste0(
  "Cases identified as any of the following categories are considered problematic and can be excluded from cSurvival analysis:"
  ,"<br>* Possible tumor/normal sample swap, cross-contamination, and/or sample purity issues;"
  ,"<br>* Case submitted is found to be a recurrence after submission;"
  ,"<br>* History of unacceptable prior treatment related to a prior/other malignancy;"
  ,"<br>* Prior malignancy;"
  ,"<br>* Item does not meet study protocol, e.g. misclassified cancer type;"
  ,"<br>* Item is noncanonical as validated by FFPE;"
  ,"<br>* Redaction."
)

# red or yellow colorscale
col_scale_no <- c(0, 0.20068666377, 0.33333333333, 0.43367666522, 0.66666666666, 1)
col_scale <- sequential_hcl(5, palette = "YlOrRd") %>% rev(.) %>% col2rgb()
col_scale <- sapply(1:ncol(col_scale), function(x) {x <- col_scale[,x]; paste0("rgb(",paste0(x, collapse = ", "),")")})
col_scale <- c("rgb(255, 255, 255)", col_scale)
col_scale <- lapply(1:length(col_scale), function(i) list(col_scale_no[i], col_scale[i]))

#TODO: ADD COLORSCALE FOR HR PLOT
#Blue and Red Orange Yellow colorscale for Hazard Ratio graph
col_scale_no <- c(0, 0.20068666377, 0.33333333333, 0.43367666522, 0.66666666666, 2)
col_scale_temp <- sequential_hcl(5, palette = "blues3") %>% col2rgb()
col_scale_hr <- sequential_hcl(5, palette = "YlOrRd") %>% rev(.) %>% col2rgb()
col_scale_hr <- append(col_scale_temp, col_scale_hr)


# col_scale <- list(list(0, "rgb(255, 255, 255)"), # 0
#                list(0.20068666377, "rgb(254,224,144)"), # 0.25 = log10(0.25) / 3
#                list(0.33333333333, "rgb(253,174,97)"), # 0.1 = log10(0.1) / 3
#                list(0.43367666522, "rgb(244,109,67)"), # 0.05 = log10(0.05) / 3
#                list(0.66666666666, "rgb(215,48,39)"), # 0.01 = log10(0.01) / 3
#                list(1, "rgb(165,0,38)") # 0.001 = log10(0.001) / 3
# )