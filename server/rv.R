rv <- reactiveValues(
  projectStatus="none", project="", max_project_n=5
  ,tcga=T,target=T,depmap=F
  
  ,variable_n_reached=0
  ,variable_n = 1
  ,verbTxtStyle1 = "box-shadow: 0 0 .2em #F5DF4D;color: red;"
  ,verbTxtStyle2 = ""
  
  ,tcga=T # check if TCGA project
  
  ,plot_type="all"
  
  ,cox_km = "cox"
  ,km_mul = "hommel" # multiple correction method
  ,median = NULL
  ,confi = T,confi_opt = "ribbon"
  ,risk_table = T,cum_table=T
  
  ,scatter_log_x=T,scatter_log_y=T,scatter_lm=T,lm_method="lm",cor_method="kendall"
  ,scatter_gender_y=F # whether to color the scatter plot by gender
)
