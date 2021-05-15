rv <- reactiveValues(
  variables_for_geo = list(), # the list to save the names and filenames of csurvival variables
  easygeo_status = FALSE, # the value to control is easyGEO is shown in our app
  show_ui = "", # "yes" upon a successful run
  
  projectStatus="none", project="", max_project_n=3, try_error=0
  ,tcga=T,target=T,depmap=F
  ,ccleStatus1="none",ccle_cancer_types="",ccle_cancer_subtypes="",depmap_gene=""
  
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
  ,censor_time_min_d=1000,censor_time_max_d=30000,censor_time_step_d=1000,censor_time_d=3000
  
  ,plot_type="all"
  
  ,cox_km = "km"
  ,km_mul = "hommel" # multiple correction method
  ,ymd="y" # default plot survival in months; "d" for days; "m" for months; "y" for years
  ,median = NULL
  ,confi = T,confi_opt = "ribbon"
  ,risk_table = T,cum_table=T
  
  ,scatter_log_x=T,scatter_log_y=T,scatter_lm=T,lm_method="lm",cor_method="kendall"
  ,scatter_gender_y=F # whether to color the scatter plot by gender
  
  # eVITTA parameters
  ,gsea_done="" # whether analyzed or not; "yes" if done
)
