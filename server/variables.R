# ------- survival analysis methods ---------
surv_methods <- c("Cox proportional-hazards (PH) model"="cox", "Kaplan-Meier (KM) log rank"="km")

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