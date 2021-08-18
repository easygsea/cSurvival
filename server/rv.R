rv <- reactiveValues(
  demo="", # yes for example run
  variables_for_geo = list(), # the list to save the names and filenames of csurvival variables
  easygeo_status = FALSE, # the value to control is easyGEO is shown in our app
  analysis_no = 0, show_ui = "", # "yes" upon a successful run
  analysis_no_hm = 0, # count of analysis for controllingheatmap plot
  flagged = "y", min_p_kc = "km",
  risk_gp = "All", risk_gpr = "All", min_gp_size = 10, risk_gps = "All",
  n_perm = 100, search_mode = "heuristic",#exhaustive
  
  projectStatus="none", project="", max_project_n=1, try_error=0
  ,tcga=T,tcgar=T,target=T,targetr=T,depmap=F,depmapr=F
  ,ccle_cancer_types="",ccle_cancer_subtypes="",depmap_gene=""
  ,depmap_path=NULL,depmap_genes=NULL,depmap_ids=NULL,depmap_ccle=NULL,cell_lines=NULL
  ,depmap_gene_appear="no"
  
  ,variable_n_reached=0
  ,variable_n = 1
  ,verbTxtStyle1 = "box-shadow: 0 0 .2em #F5DF4D;color: red;"
  ,verbTxtStyle2 = ""
  
  ,tcga=T # check if TCGA project
  ,tcga_stype="OS" # if TCGA, type of survival analysis
  ,plot_stype="OS" # label to display
  ,tcga_code="" # clinical data quality info for the selected TCGA project
  ,tcga_msg="" # warning message if NA annotation for TCGA clinical outcome codes; "" means no NA
  ,censor_time_ymd="y" # default censor time unit: years
  ,censor_time=10 # default censor time 10 yrs
  ,censor_time_min=1,censor_time_max=100,censor_time_step=1
  ,censor_time_min_y=1,censor_time_max_y=100,censor_time_step_y=1,censor_time_y=10
  ,censor_time_min_m=20,censor_time_max_m=2000,censor_time_step_m=20,censor_time_m=200
  ,censor_time_min_d=400,censor_time_max_d=30000,censor_time_step_d=200,censor_time_d=3000
  
  ,plot_type="all"
  
  ,cox_km = "km",cox_kmr="km"
  ,km_mul = "hommel" # multiple correction method
  ,km_mul_padj = "padj"
  ,km_mul_dp = "hommel"
  ,km_mul_dp_padj = "padj"
  ,ymd="y" # default plot survival in months; "d" for days; "m" for months; "y" for years
  ,ymd_int=1,ymd_int_range=c(1,2,3)
  ,ymd_int_y=1,ymd_int_range_y=c(1,2,3)
  ,ymd_int_m=20,ymd_int_range_m=c(1,2,3,5,10,20,40)
  ,ymd_int_d=1000,ymd_int_range_d=c(30,50,100,200,300,1000)
  ,median = NULL
  ,confi = T,confi_opt = "ribbon"
  ,risk_table = T,cum_table=T
  
  ,dens_fill=T, dens_mean=T
  
  ,palette="jco"
  
  ,scatter_log_x=T,scatter_log_y=T,scatter_lm=T,lm_method="lm",cor_method="kendall",sm_conf=T
  ,scatter_gender_y=F # whether to color the scatter plot by gender
  
  # eVITTA parameters
  ,gsea_done="" # whether analyzed or not; "yes" if done
  
  # source data options
  ,project_tcga="TCGA-ACC"
  
  # hr & p value VS quantile graph
  ,quantile_graph=NULL
  
  # annotate cells
  ,annot_cells_y="" # yes
  # annotate data points
  ,annot_data_points_y=""
  
  # # permutations
  # ,nitem = 100
  
  ,violin_k=1.5 # no of s.d. to plot on violin
  ,violin_log_y=T
  ,violin_trim=T
)
