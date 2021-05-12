rv <- reactiveValues(
  variables_for_geo = list(), # the list to save the names and filenames of csurvival variables
  easygeo_status = FALSE, # the value to control is easyGEO is shown in our app
  show_ui = "", # "yes" upon a successful run
  
  projectStatus="none", project="", max_project_n=3, try_error=0
  ,tcga=T,target=T,depmap=F
  
  ,variable_n_reached=0
  ,variable_n = 1
  ,verbTxtStyle1 = "box-shadow: 0 0 .2em #F5DF4D;color: red;"
  ,verbTxtStyle2 = ""
  
  ,tcga=T # check if TCGA project
  ,tcga_stype="OS" # if TCGA, type of survival analysis
  ,plot_stype="OS" # label to display
  ,censor_time="10" # default censor time 10 yrs

  ,plot_type="all"
  
  ,cox_km = "cox"
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
