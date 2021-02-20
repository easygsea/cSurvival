# dynamic variables need initiatialization and updates according to inputs' changes
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
    ,gs_lib_genes_id <- paste0("gs_lgs_",x)
    ,gs_gene_id <- paste0("gs_lg_",x)
    # manual gene input
    ,gs_manual_id <- paste0("gs_m_",x)
    ,gs_genes_id <- paste0("gs_mg_",x)
    ,lower_id <- paste0("lower_",x)
    ,higher_id <- paste0("upper_",x)
    ,step_id <- paste0("step_",x)
  )
}